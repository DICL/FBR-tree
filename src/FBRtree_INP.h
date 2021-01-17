#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "sys/time.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <pthread.h>
#include <queue>	// using queue C++ style
#include <stack>	// using queue C++ style
#include <cstring>
#include "forPM.h"
#include "index.h"

#include <shared_mutex>
#include <mutex>

#define BIG_NUM (FLT_MAX/4.0)

#define Undefined(x) ((x)->boundary[0] > (x)->boundary[NUMDIMS])
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#define cacheLineSize 64
extern int SEARCH;
extern int IpthreadNum;
extern uint64_t flipCount;
extern uint64_t k_tmp_flipCount;
extern uint64_t k_split_flipCount;
extern uint64_t clflushCnt;
extern uint64_t k_tmp_clflushCnt;
extern uint64_t k_split_clflushCnt;


#ifdef BREAKDOWN
double traversal_time = 0.0;
double write_time = 0.0;
struct timeval tr1,tr2;

extern double k_insert_time;
#ifdef SPLIT
double split_write_time = 0.0;
extern double k_split_insert_time;
int split_cnt = 0;
#endif
#endif
int isSplit = 0;

void getClf(char* data, int len){
	volatile char *ptr = (char *)((unsigned long)data &~(cacheLineSize-1));
	for(; ptr<data+len; ptr+=cacheLineSize){
		k_tmp_clflushCnt++;
	}
}

// using queue C++ style
using namespace std;
queue<struct Node*> nodeQueue;

int Compare(struct Rect *o, struct Rect *n);
void RTreePrint(TOID(Node) n);

TOID(Node) RTreeNewNode(PMEMobjpool*);
struct Rect RTreeNodeCover(TOID(Node) n);
int RTreePickBranch(struct Rect *r, TOID(Node) n);
void RTreeInitRect(struct Rect *r);
float RTreeRectVolume(struct Rect *r);
struct Rect RTreeCombineRect(struct Rect *r, struct Rect *rr);
inline int RTreeOverlap(struct Rect *r, struct Rect *s);
int RTreeContained(struct Rect *r, struct Rect *s);
void RTreeSplitNode(PMEMobjpool *pop, TOID(Node) n, struct Branch *b, TOID(Node) *nn, TOID(Node) p);

double checkFreeSpace(TOID(Node) n)
{
	uint64_t freeSpace=0;
	uint64_t nfreeSpace=0;
	int nodeCount=0;
	int nodeSeq = 0;
	nodeQueue.push((Node*)n.oid.off);
	nodeSeq++;

	while(!nodeQueue.empty()){
		TOID(Node) t = n;
		t.oid.off = (uint64_t)nodeQueue.front();
		nodeQueue.pop();
		nodeCount++;
		if(!D_RW(t)->meta.IsLeaf()){ // This is an internal node in the tree
			for(int i=0; i<NODECARD; i++) {	
				struct Branch b = D_RW(t)->branch[i];
				if(D_RW(t)->meta.Bit(i)){ 
					nodeQueue.push(b.child);
					nodeSeq++;
					nfreeSpace++;
				}
				else
					freeSpace++;
			}
		} else {	
			for(int i=0; i<NODECARD; i++) {
				if(!D_RW(t)->meta.Bit(i))
					freeSpace++;
				else 
					nfreeSpace++;
			}
		}
	}
	return freeSpace*sizeof(Branch) / nodeCount*NODECARD*sizeof(Branch);
}

void RTreePrint(TOID(Node) n)
{
	int nodeSeq = 0;
	nodeQueue.push((Node*)n.oid.off);
	nodeSeq++;

	while(!nodeQueue.empty()) {
		TOID(Node) t = n;
		t.oid.off = (uint64_t)nodeQueue.front();
		nodeQueue.pop();

    if(!D_RW(t)->meta.IsLeaf()) { // This is an internal node in the tree
			D_RW(t)->meta.Print();
      printf("------Not leaf : [%p]-------\n", t.oid.off);
      for(int i=0; i<NODECARD; i++) {
				struct Branch b = D_RO(t)->branch[i];
		
				if(D_RW(t)->meta.Bit(i)) {
          printf("%lf %lf %lf %lf %f %f\n", b.rect.boundary[0], b.rect.boundary[1], 
						b.rect.boundary[2], b.rect.boundary[3], b.rect.boundary[4], b.rect.boundary[5]); 
					nodeQueue.push(b.child);
					nodeSeq++;
				}
			}
			printf("\n");
		} else {	
			printf("------ Leaf : [%p]-------\n", t.oid.off);
			D_RO(t)->meta.Print();
			for(int i=0; i<NODECARD; i++) {
				struct Branch b = D_RO(t)->branch[i];
				if(D_RW(t)->meta.Bit(i))
					printf("%lf %lf %lf %lf %f %f\n", b.rect.boundary[0], b.rect.boundary[1], 
						b.rect.boundary[2], b.rect.boundary[3], b.rect.boundary[4], b.rect.boundary[5]); 
			}
			printf("\n");
		}
	}
}

// Initialize one branch cell in a node.
static void RTreeInitBranch(struct Branch *b)
{
	register int i;
	RTreeInitRect(&(b->rect));
	b->child = NULL;
}

// Initialize a Node structure.
void RTreeInitNode(TOID(Node) n)
{
	register int i;
	D_RW(n)->mutex_ = new std::shared_mutex();
	D_RW(n)->meta.Init();
	D_RW(n)->meta.Leaf();
	for(int i=0; i<NODECARD; i++)
		RTreeInitBranch(&(D_RW(n)->branch[i])); 
}

// Initialize the new r-tree node
TOID(Node) RTreeNewNode(PMEMobjpool *pop)
{
	TOID(Node) n = OID_NULL;
	POBJ_ALLOC(pop, &n, Node, sizeof(Node), NULL, NULL);

	RTreeInitNode(n);
	return n;
}

// Initialize the r-tree new index
TOID(Node) RTreeNewIndex(PMEMobjpool *pop, TOID(Node) root)
{
	RTreeInitNode(root);
	/* leaf */
	return root;
} 

extern TOID(Node) total_root;

thread_local TOID(splitLog) thread_log;

// R-tree search in CPU
uint64_t hostRTreeSearch(TOID(Node) n, struct Rect *r)
{
  uint64_t hitCount = 0;
  int i;
  uint64_t ret = 0;
  uint8_t version;
#ifdef SHARED
  D_RW(n)->mutex_->lock_shared();
#endif

  //lock-free 
  if (!D_RW(n)->meta.IsLeaf()) { // this is an internal node in the tree 
    version = D_RO(n)->meta.Version();
    std::queue<Node*> childqueue;

    while(1) {
      for (i=0; i<NODECARD; i++) {
        //get all overlap branches
        if (RTreeOverlap(r,&D_RW(n)->branch[i].rect) && D_RW(n)->meta.Bit(i)) {
          TOID(Node) child;
          child.oid.off = (uint64_t)D_RW(n)->branch[i].child;	
          childqueue.push((Node*)child.oid.off);
        }
      }

      if(version == D_RW(n)->meta.Version() || D_RW(n)->meta.Version() == 0) {
        while(!childqueue.empty()) {
          TOID(Node) nextNode = n;
          nextNode.oid.off = (uint64_t)childqueue.front();
          ret = hostRTreeSearch(nextNode, r);

          if (ret == -1) {
            std::queue<Node*> newchild;
            std::swap(childqueue, newchild);
            version = D_RW(n)->meta.Version();
            break;
          }

          hitCount += ret;
          childqueue.pop();
        }

        if (ret == -1)
          continue;

        break;
      }

      else { // child split and added branch
        // bresion changed, it means new branch added
        std::queue<Node*> newchild;
        std::swap(childqueue, newchild);
        version = D_RW(n)->meta.Version();            
      }
    }
  }

  else {
    // this is a lead node
    version = D_RW(n)->meta.Version();
    uint64_t leaf_cnt = 0;
    for (i=0; i<NODECARD; i++) {
      if (!D_RW(n)->meta.Bit(i))
        continue;
      if (RTreeOverlap(r, &D_RW(n)->branch[i].rect))
        leaf_cnt++;
    }
    if (version == D_RW(n)->meta.Version()) {
#ifdef SHARED
      D_RW(n)->mutex_->unlock_shared();
#endif  
      return leaf_cnt;
    } else {
#ifdef SHARED
      D_RW(n)->mutex_->unlock_shared();
#endif
      return -1;
    }
  }
#ifdef SHARED
  D_RW(n)->mutex_->unlock_shared();
#endif
  return hitCount;
}

// Inserts a new data rectangle into the index structure.
// Recursively descends tree, propagates splits back up.
// Returns 0 if node was not split.  Old node updated.
// If node was split, returns 1 and sets the pointer pointed to by
// new_node to point to the new node.  Old node updated to become one of two.
// The level argument specifies the number of steps up from the leaf
// level to insert; e.g. a data rectangle goes in at level = 0.
//
thread_local int LogCur = 0;

static int RTreeInsertRect2(PMEMobjpool *pop, Rect *r, uint64_t dataID, TOID(Node) n, TOID(Node) *new_node, TOID(Node) p, TOID(splitLog) sl)
{
  int i;
  struct Branch b;

  //n_ is child holder
  TOID(Node) n_;
  TOID(Node) n2;
  TOID(Node) parent = p;
  int rt = 0;
  bool unlocked = false;
  // Still above level for insertion, go down tree recursively
  //

  if (!D_RW(n)->meta.IsLeaf()) {
    i = RTreePickBranch(r, n); // i is the index of split branc
    parent = n;

#ifdef BREAKDOWN
    gettimeofday(&tr2,0); // start the stopwatch
    traversal_time += (tr2.tv_sec-tr1.tv_sec)*1000000 + (tr2.tv_usec - tr1.tv_usec);
#endif

    D_RW(n)->branch[i].rect = RTreeCombineRect(r, &D_RW(n)->branch[i].rect);
    pmemobj_persist(pop, (char*)&D_RO(n)->branch[i], sizeof(Branch));
    getClf((char*)&D_RO(n)->branch[i], sizeof(Branch));
    k_tmp_flipCount += sizeof(struct Branch);

#ifdef BREAKDOWN
    gettimeofday(&tr1,0); // start the stopwatch
#endif

    TOID(Node) childn = n;
    childn.oid.off = (uint64_t)D_RO(n)->branch[i].child;
    D_RW(childn)->mutex_->lock();

    if(!D_RO(childn)->meta.IsFull()) {
      unlocked = true;
      D_RW(n)->mutex_->unlock();
    }

    //Store child in n_
    n_ = n;
    n_.oid.off = (uint64_t)D_RW(n)->branch[i].child;

    if (!RTreeInsertRect2(pop, r, dataID, n_, &n2, n, OID_NULL)) { 
      // child was not split
      if(!unlocked) {
        D_RW(n)->mutex_->unlock();
      }
      return 0;
    }
    else {
      // child was split
      isSplit = 1;

      if(unlocked) {
        D_RW(n)->mutex_->lock();
      }
      b.child = (Node*) n2.oid.off;
      b.rect = RTreeNodeCover(n2);

#ifdef FULLLOG
      POBJ_FREE(&thread_log);
#endif
      rt = RTreeAddBranch(pop, &b, n, new_node, p, OID_NULL);

      // child->VersionIncr();
      D_RW(n_)->meta.VersionIncr();
      pmemobj_persist(pop, (char*)D_RO(n_), 1);
      getClf((char*)D_RO(n_), 1);
      k_tmp_flipCount += 1;

      n_.oid.off = (uint64_t)D_RW(n)->branch[i].child;
      D_RW(n)->branch[i].rect = RTreeNodeCover(n_);
      pmemobj_persist(pop, (char*)&D_RW(n)->branch[i], sizeof(struct Branch));
      getClf((char*)&D_RW(n)->branch[i], sizeof(struct Branch));
      k_tmp_flipCount += sizeof(Branch);

      if(!(rt == 1 && n.oid.off == p.oid.off)) {
        D_RW(n)->mutex_->unlock();
      }
    }
    return rt;
  } else {
    // Have reached level for insertion. Add rect, split if necessary
#ifdef BREAKDOWN
    gettimeofday(&tr2,0); // start the stopwatch
    traversal_time += (tr2.tv_sec-tr1.tv_sec)*1000000 + (tr2.tv_usec - tr1.tv_usec);
#endif

    b.rect = *r;
    b.child = (Node *) dataID;
    /* child field of leaves contains dataID of data record */
    rt = RTreeAddBranch(pop, &b, n, new_node, p, OID_NULL);
    if(!(rt == 1 && n.oid.off == p.oid.off)) {
      D_RW(n)->mutex_->unlock();
    } 
    return rt;
  }
}

// Insert a data rectangle into an index structure.
// RTreeInsertRect provides for splitting the root;
// returns 1 if root was split, 0 if it was not.
// The level argument specifies the number of steps up from the leaf
// level to insert; e.g. a data rectangle goes in at level = 0.
// RTreeInsertRect2 does the recursion.
//
int RTreeInsertRect(PMEMobjpool *pop, struct Rect *R, uint64_t Did, TOID(Node) Root, TOID(splitLog) sl)
{
	register struct Rect *r = R;
	register uint64_t dataID = Did;
	TOID(Node) root = Root;
	TOID(Node) newroot;

	TOID(Node) newnode = OID_NULL;
	struct Branch b;
	int result;

#ifdef BREAKDOWN
	gettimeofday(&tr1,0); // start the stopwatch
#endif

	TOID(Node) my_root = root;
	while(1) {
		D_RW(my_root)->mutex_->lock();
		if(my_root.oid.off == root.oid.off) {
			break;
		}
		D_RW(my_root)->mutex_->unlock();
		my_root = root;
	}

	if (RTreeInsertRect2(pop, r, dataID, root, &newnode, root, OID_NULL)) {
    /* root split */
    isSplit = 1;
    //make new root
		newroot = RTreeNewNode(pop);  /* grow a new root, & tree taller */

#ifdef FULLLOG
    POBJ_FREE(&thread_log);
    POBJ_ALLOC(pop, &thread_log, splitLog, sizeof(splitLog), NULL, NULL);
    
    D_RW(thread_log)->parentPoint = (Node*)newroot.oid.off;
    D_RW(thread_log)->childPoint = (Node*)root.oid.off;
    D_RW(thread_log)->siblingPoint = (Node*)newnode.oid.off;

    memcpy((void*)&D_RW(thread_log)->parent, D_RW(newroot), sizeof(Node));
    memcpy((void*)&D_RW(thread_log)->child, D_RW(root), sizeof(Node));
    memcpy((void*)&D_RW(thread_log)->sibling, D_RW(newnode), sizeof(Node));

    pmemobj_persist(pop, (char*)D_RW(thread_log), sizeof(struct splitLog));
    getClf((char*)D_RW(thread_log), sizeof(struct splitLog));
    k_tmp_flipCount += sizeof(splitLog);
#endif

		//make new branch which point original root
		D_RW(newroot)->meta.Iter();
		b.rect = RTreeNodeCover(root);
		b.child = (Node*)root.oid.off;
		RTreeAddBranch(pop, &b, newroot, 0, OID_NULL, OID_NULL);

		b.rect = RTreeNodeCover(newnode);
		b.child = (Node*)newnode.oid.off;
		RTreeAddBranch(pop, &b, newroot, 0, OID_NULL, OID_NULL);
		D_RW(newroot)->meta.VersionIncr();

		pmemobj_persist(pop, (char *)D_RO(newroot), META + 2 * sizeof(Branch));
		getClf((char *)D_RO(newroot), META + 2 * sizeof(Branch));
		k_tmp_flipCount += 8 + 2*sizeof(Branch);
		
		TOID(Node) temp = root;
		D_RW(temp)->meta.VersionIncr();
		pmemobj_persist(pop, (char*)D_RO(temp),1);
		k_tmp_flipCount += 1;
		getClf((char*)D_RO(temp),1);
		
		total_root.oid.off = newroot.oid.off;

		{   
			D_RW(newroot)->mutex_->lock();
			{
				root = newroot;
				D_RW(temp)->mutex_->unlock();
			}
			D_RW(newroot)->mutex_->unlock();
		}
    // free rootLog
		POBJ_FREE(&thread_log);
    result = 1;
	} else {
		result = 0;
	}
	return result;
}

// Find the smallest rectangle that includes all rectangles in
// branches of a node.
//
struct Rect RTreeNodeCover(TOID(Node) n)
{
	int i, first_time=1;
	struct Rect r;

	RTreeInitRect(&r);
	for (i = 0; i < NODECARD; i++) {
		if (D_RW(n)->meta.Bit(i)) {
			if (first_time) {
				r = D_RO(n)->branch[i].rect;
				first_time = 0;
			}
			else
				r = RTreeCombineRect(&r, &(D_RW(n)->branch[i].rect));
		}
	}
	return r;
}

// Pick a branch.  Pick the one that will need the smallest increase
// in area to accomodate the new rectangle.  This will result in the
// least total area for the covering rectangles in the current node.
// In case of a tie, pick the one which was smaller before, to get
// the best resolution when 	int SEARCH = atoi(args[3]);searching.
//
int RTreePickBranch(struct Rect *r, TOID(Node) n)
{
	struct Rect *rr;
	int i = 0, first_time=1;

	float increase, bestIncr=(float)-1, area, bestArea;
	int best = -1;
	struct Rect tmp_rect;

	for (i=0; i<NODECARD; i++) {
		if(!D_RW(n)->meta.Bit(i)) 
            continue;

		rr = &D_RW(n)->branch[i].rect;
		area = RTreeRectVolume(rr);
		tmp_rect = RTreeCombineRect(r, rr);
		increase = RTreeRectVolume(&tmp_rect) - area;

    if(increase == 0) {
			return i;
		}
	
    if (increase < bestIncr || first_time) {
			best = i;
			bestArea = area;
			bestIncr = increase;
			first_time = 0;
		} else if (increase == bestIncr && area < bestArea) {
			best = i;
			bestArea = area;
			bestIncr = increase;
		}
	}
	return best;
}


// Add a branch to a node.  Split the node if necessary.
// Returns null if node not split.  Old node updated.
// Returns a pointer to a new node if node splits, sets *new_node to address of new node.
// Old node updated, becomes one of two.
//
int RTreeAddBranch(PMEMobjpool *pop, struct Branch *B, TOID(Node) N, TOID(Node) *New_node, TOID(Node) PN, TOID(splitLog) sl)
{
	register struct Branch *b = B;
	TOID(Node) n = N;
	TOID(Node) *new_node = New_node;
	TOID(Node) p = PN;
	register int i, j;
	
	if (!D_RO(n)->meta.IsFull()) {
    /* split won't be necessary */
#ifdef FULLLOG
    TOID(Node) nodeLog;
    POBJ_ALLOC(pop, &nodeLog, Node, sizeof(Node), NULL, NULL);
    pmemobj_memcpy_persist(pop, (char*)D_RW(nodeLog), D_RW(n), sizeof(Node));
    getClf((char*)D_RW(nodeLog), sizeof(Node));
    k_tmp_flipCount += sizeof(Node);
#endif

    for (i = 0; i < NODECARD; i++) {  /* find empty branch */
			if (!D_RW(n)->meta.Bit(i)) {
				D_RW(n)->branch[i] = *b;
				D_RW(n)->meta.Set(i);
        k_tmp_flipCount+=1;
        pmemobj_persist(pop, (char*)&D_RW(n)->branch[i], sizeof(struct Branch));
        k_tmp_flipCount += sizeof(Branch);
				getClf((char*)&D_RW(n)->branch[i], sizeof(struct Branch));
        D_RW(n)->meta.VersionIncr();
        pmemobj_persist(pop, (char*)D_RW(n), META);
        k_tmp_flipCount += META;
				getClf((char*)D_RW(n), META);
				break;
			}
		}
#ifdef FULLLOG
    POBJ_FREE(&nodeLog);
#endif
		return 0;
	} else {
    *new_node = RTreeNewNode(pop);
    TOID(Node) nn = *new_node;

#ifdef FULLLOG
    POBJ_ALLOC(pop, &thread_log, splitLog, sizeof(splitLog), NULL, NULL);
    D_RW(thread_log)->parentPoint = (Node*)p.oid.off;
    D_RW(thread_log)->childPoint = (Node*)n.oid.off;
    D_RW(thread_log)->siblingPoint = (Node*)nn.oid.off;

    memcpy((void*)&D_RW(thread_log)->parent, D_RW(p), sizeof(Node));
    memcpy((void*)&D_RW(thread_log)->child, D_RW(n), sizeof(Node));
    memcpy((void*)&D_RW(thread_log)->sibling, D_RW(nn), sizeof(Node));

    pmemobj_persist(pop, (char*)D_RW(thread_log), sizeof(struct splitLog));
    getClf((char*)D_RW(thread_log), sizeof(struct splitLog));
    k_tmp_flipCount += sizeof(splitLog);
#endif

    RTreeSplitNode(pop, n, b, new_node, p);

    return 1;
  }
}

/*-----------------------------------------------------------------------------
| Initialize a rectangle to have all 0 coordinates.
-----------------------------------------------------------------------------*/
void RTreeInitRect(struct Rect *r)
{
	int i;
	for (i=0; i<NUMSIDES; i++)
		r->boundary[i] = (float)0;
}

/*-----------------------------------------------------------------------------
| Calculate the n-dimensional volume of a rectangle
-----------------------------------------------------------------------------*/
float RTreeRectVolume(struct Rect *r)
{
	int i;
	float volume = (float)1;

	if (Undefined(r))
		return (float)0;

	for(i=0; i<NUMDIMS; i++)
		volume *= r->boundary[i+NUMDIMS] - r->boundary[i];

	return volume;
}

/*-----------------------------------------------------------------------------
| Combine two rectangles, make one that includes both.
----------------------------------------false-------------------------------------*/
struct Rect RTreeCombineRect(struct Rect *r, struct Rect *rr)
{
	int i, j;
	struct Rect new_rect;

	if (Undefined(r))
		return *rr;

	if (Undefined(rr))
		return *r;

	for (i = 0; i < NUMDIMS; i++) {
		new_rect.boundary[i] = MIN(r->boundary[i], rr->boundary[i]);
		j = i + NUMDIMS;
		new_rect.boundary[j] = MAX(r->boundary[j], rr->boundary[j]);
	}
	return new_rect;
}

int Compare(struct Rect *o, struct Rect *n)
{
	int num = NUMDIMS*2;
	for (int i = 0; i <NUMDIMS*2; i++) {
		if(o->boundary[i] == n->boundary[i]) {
			num--;
        }
	}
	return num;	
}

/*-----------------------------------------------------------------------------
| Decide whether two rectangles overlap.
-----------------------------------------------------------------------------*/
inline int RTreeOverlap(struct Rect *r, struct Rect *s)
{
	int i, j;

	for (i=0; i<NUMDIMS; i++) {
		j = i + NUMDIMS;  
		if (r->boundary[i] > s->boundary[j] ||
				s->boundary[i] > r->boundary[j]) {
			return FALSE;
		}
	}
	return TRUE;
}

/*-----------------------------------------------------------------------------
| Decide whether rectangle r is contained in rectangle s.
-----------------------------------------------------------------------------*/
int RTreeContained(struct Rect *r, struct Rect *s)
{
	int i, j, result;

	// undefined rect is contained in any other
	//
	if (Undefined(r))
		return TRUE;

	// no rect (except an undefined one) is contained in an undef rect
	//
	if (Undefined(s))
		return FALSE;

	result = TRUE;
	for (i = 0; i < NUMDIMS; i++) {
		j = i + NUMDIMS;  /* index for high sides */
		result = result
			&& r->boundary[i] >= s->boundary[i]
			&& r->boundary[j] <= s->boundary[j];
	}
	return result;
}

/*-----------------------------------------------------------------------------
| Load branch buffer with branches from full node plus the extra branch.
-----------------------------------------------------------------------------*/
static void RTreeGetBranches(struct forSplit *fs, TOID(Node) n, struct Branch *b)
{
	int i;

	/* load the branch buffer */
	for (i=0; i<NODECARD; i++) {
		fs->BranchBuf[i] = D_RO(n)->branch[i];
	}
	fs->BranchBuf[NODECARD] = *b;
	fs->BranchCount = NODECARD + 1;

	/* calculate rect containing all in the set */
	fs->CoverSplit = fs->BranchBuf[0].rect;
	for (i=1; i<NODECARD+1; i++) {
		fs->CoverSplit = RTreeCombineRect(&fs->CoverSplit, &fs->BranchBuf[i].rect);
	}	
}


/*-----------------------------------------------------------------------------
| Put a branch in one of the groups.
-----------------------------------------------------------------------------*/
static void RTreeClassify(struct forSplit* fs, int i, int group, struct PartitionVars *p)
{
	p->partition[i] = group;
	p->taken[i] = TRUE;

	if (p->count[group] == 0)
		p->cover[group] = fs->BranchBuf[i].rect;
	else
		p->cover[group] = RTreeCombineRect(&fs->BranchBuf[i].rect,
				&p->cover[group]);
	p->area[group] = RTreeRectVolume(&p->cover[group]);
	p->count[group]++;
}

/*-----------------------------------------------------------------------------
| Initialize a PartitionVars structure.
-----------------------------------------------------------------------------*/
static void RTreePickSeeds(struct forSplit* fs, struct PartitionVars *P)
{
	register struct PartitionVars *p = P;
	register int i, dim, high;
	register struct Rect *r, *rlow, *rhigh;
	register float w, separation, bestSep;
	RectReal width[NUMDIMS];
	int leastUpper[NUMDIMS], greatestLower[NUMDIMS];
	int seed0, seed1;

	for (dim=0; dim<NUMDIMS; dim++) {
		high = dim + NUMDIMS;

		/* find the rectangles farthest out in each direction
		 * along this dimens */
		greatestLower[dim] = leastUpper[dim] = 0;
		for (i=1; i<NODECARD+1; i++) {
			r = &fs->BranchBuf[i].rect;
			if (r->boundary[dim] >
					fs->BranchBuf[greatestLower[dim]].rect.boundary[dim]) {
				greatestLower[dim] = i;
			}
			if (r->boundary[high] <
					fs->BranchBuf[leastUpper[dim]].rect.boundary[high]) {
				leastUpper[dim] = i;
			}
		}

		/* find width of the whole collection along this dimension */
		width[dim] = fs->CoverSplit.boundary[high] -
			fs->CoverSplit.boundary[dim];
	}

	/* pick the best separation dimension and the two seed rects */
	for (dim=0; dim<NUMDIMS; dim++)	{
		high = dim + NUMDIMS;

		/* divisor for normalizing by width */
		if (width[dim] == 0)
			w = (RectReal)1;
		else
			w = width[dim];

		rlow = &fs->BranchBuf[leastUpper[dim]].rect;
		rhigh = &fs->BranchBuf[greatestLower[dim]].rect;
		if (dim == 0) {
			seed0 = leastUpper[0];
			seed1 = greatestLower[0];
			separation = bestSep =
				(rhigh->boundary[0] -
				 rlow->boundary[NUMDIMS]) / w;
		} else {
			separation =
				(rhigh->boundary[dim] -
				 rlow->boundary[dim+NUMDIMS]) / w;
			if (separation > bestSep) {
				seed0 = leastUpper[dim];
				seed1 = greatestLower[dim];
				bestSep = separation;
			}
		}
	}

	if (seed0 != seed1) {
		RTreeClassify(fs, seed0, 0, p);
		RTreeClassify(fs, seed1, 1, p);
	}
}

/*-----------------------------------------------------------------------------
| Put each rect that is not already in a group into a group.
| Process one rect at a time, using the following hierarchy of criteria.
| In case of a tie, go to the next test.
| 1) If one group already has the max number of elements that will allow
| the minimum fill for the other group, put r in the other.
| 2) Put r in the group whose cover will expand less.  This automatically
| takes care of the case where one group cover contains r.
| 3) Put r in the group whose cover will be smaller.  This takes care of the
| case where r is contained in both covers.
| 4) Put r in the group with fewer elements.
| 5) Put in group 1 (arbitrary).
|
| Also update the covers for both groups.
-----------------------------------------------------------------------------*/
static void RTreePigeonhole(struct forSplit* fs, struct PartitionVars *P)
{
	register struct PartitionVars *p = P;
	struct Rect newCover[2];
	register int i, group;
	RectReal newArea[2], increase[2];

	for (i=0; i<NODECARD+1; i++) {
		if (!p->taken[i]) {
			/* if one group too full, put rect in the other */
			if (p->count[0] >= p->total - p->minfill) {
				RTreeClassify(fs, i, 1, p);
				continue;
			} else if (p->count[1] >= p->total - p->minfill) {
				RTreeClassify(fs, i, 0, p);
				continue;
			}

			/* find areas of the two groups' old and new covers */
			for (group=0; group<2; group++) {
				if (p->count[group]>0)
					newCover[group] = RTreeCombineRect(
							&fs->BranchBuf[i].rect,
							&p->cover[group]);
				else
					newCover[group] = fs->BranchBuf[i].rect;
				newArea[group] = RTreeRectVolume(
						&newCover[group]);
				increase[group] = newArea[group]-p->area[group];
			}

			/* put rect in group whose cover will expand less */
			if (increase[0] < increase[1])
				RTreeClassify(fs, i, 0, p);
			else if (increase[1] < increase[0])
				RTreeClassify(fs, i, 1, p);

			/* put rect in group that will have a smaller cover */
			else if (p->area[0] < p->area[1])
				RTreeClassify(fs, i, 0, p);
			else if (p->area[1] < p->area[0])
				RTreeClassify(fs, i, 1, p);

			/* put rect in group with fewer elements */
			else if (p->count[0] < p->count[1])
				RTreeClassify(fs, i, 0, p);
			else
				RTreeClassify(fs, i, 1, p);
		}
	}
}

/*-----------------------------------------------------------------------------
| Pick two rects from set to be the first elements of the two groups.
| Pick the two that are separated most along any dimension, or overlap least.
| Distance for separation or overlap is measured modulo the width of the
| space covered by the entire set along that dimension.
-----------------------------------------------------------------------------*/
static void RTreeInitPVars(struct PartitionVars *P, int maxrects, int minfill)
{
	register struct PartitionVars *p = P;
	register int i;

	p->count[0] = p->count[1] = 0;
	p->total = maxrects;
	p->minfill = minfill;
	for (i=0; i<maxrects; i++) {
		p->taken[i] = FALSE;
		p->partition[i] = -1;
	}
}

/*-----------------------------------------------------------------------------
| Method 0 for finding a partition:
| First find two seeds, one for each group, well separated.
| Then put other rects in whichever group will be smallest after addition.
-----------------------------------------------------------------------------*/
static void RTreeMethodZero(struct forSplit* fs, struct PartitionVars *p, int minfill)
{
	RTreeInitPVars(p, fs->BranchCount, minfill);
	RTreePickSeeds(fs, p);
	RTreePigeonhole(fs, p);
}

/*-----------------------------------------------------------------------------
| Copy branches from the buffer into two nodes according to the partition.
-----------------------------------------------------------------------------*/
static void RTreeLoadNodes(PMEMobjpool *pop, struct forSplit* fs, TOID(Node) N, TOID(Node) Q,
    struct PartitionVars *P)
{
  TOID(Node) n = N, q = Q;
  register struct PartitionVars *p = P;
  register int i;


  int newA = p->partition[NODECARD];

  if(D_RW(n)->meta.IsLeaf())
    D_RW(q)->meta.Leaf();
  else
    D_RW(q)->meta.Iter();
  D_RW(n)->meta.VersionReset();
  D_RW(q)->meta.VersionIncr();
  int count = 0;
  for (i=0; i<NODECARD+1; i++) {
    if (p->partition[i] == newA) {
      D_RW(q)->branch[count] = fs->BranchBuf[i];;
      D_RW(q)->meta.Set(count);
      if(i < NODECARD){
        D_RW(n)->meta.Reset(i);
        flipCount++;
      }
      flipCount++;
      count++;
    }
  }
  pmemobj_persist(pop, (char*)D_RO(q), META + sizeof(Branch) * count);
  getClf((char*)D_RO(q), META + sizeof(Branch) * count);
  k_tmp_flipCount += META + sizeof(Branch) * count;
  pmemobj_persist(pop, (char*)D_RO(n), META); 
  getClf((char*)D_RO(n), META);
  k_tmp_flipCount += META;
}

/*-----------------------------------------------------------------------------
| Split a node.
| Divides the nodes branches and the extra one between two nodes.
| Old node is one of the new ones, and one really new one is created.
-----------------------------------------------------------------------------*/
void RTreeSplitNode(PMEMobjpool *pop, TOID(Node) n, struct Branch *b, TOID(Node) *nn, TOID(Node) pn)
{
	register struct forSplit fs;
  register struct PartitionVars *p;
	
	/* load all the branches into a buffer, initialize old node */
	RTreeGetBranches(&fs, n, b);
	/* find partition */
	p = fs.Partitions;

	/* Note: can't use MINFILL(n) below since n was cleared by GetBranches() */
	RTreeMethodZero(&fs, p, NODECARD/2 );


	/* record how good the split was for statistics */

	/* put branches from buffer in 2 nodes according to chosen partition */
	RTreeLoadNodes(pop, &fs, n, *nn, p);
}
