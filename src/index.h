#ifndef _INDEX_
#define _INDEX_

#include <stdbool.h>
#include <libpmemobj.h>
#include <stdint.h>
#include <shared_mutex>
#define NUMDIMS	3	/* number of dimensions */
#define NDEBUG

typedef float RectReal;
typedef bool DeadFlag;

/*-----------------------------------------------------------------------------
| Global definitions.
-----------------------------------------------------------------------------*/

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

#define NUMSIDES 2*NUMDIMS // Size = 6

typedef struct Rect Rect;
typedef struct Node Node;
typedef struct Branch Branch;
typedef struct splitLog splitLog;
typedef struct cowLog cowLog;

POBJ_LAYOUT_BEGIN(rtree);
POBJ_LAYOUT_ROOT(rtree, Node);
POBJ_LAYOUT_TOID(rtree, Rect);
POBJ_LAYOUT_TOID(rtree, Branch);
POBJ_LAYOUT_TOID(rtree, splitLog);
POBJ_LAYOUT_TOID(rtree, cowLog);
POBJ_LAYOUT_END(rtree);

struct Rect // Float(4byte) * 6
{
	RectReal boundary[NUMSIDES]; /* xmin,ymin,...,xmax,ymax,... */
};

struct Node;

#ifndef MULTIMETA
#define META 8
#define DUMMY 24
#define MAXCARD 55
#define NODECARD 55
#else
#ifdef NS2
#define META 16
#define DUMMY 16
#define MAXCARD 119
#define NODECARD 119
#elif NS3
#define META 24
#define DUMMY 8
#define MAXCARD 183
#define NODECARD 183
#elif NS4
#define META 32
#define DUMMY 0
#define MAXCARD 247
#define NODECARD 247
#endif
#endif
#include "bitmap.h" 
//
struct Branch // 24 + 8
{
	Rect rect; // Float(4byte) * 6
	Node* child; // Pointer =  8
};

/* max branching factor of a node */
struct Node
{
    MetaData<META> meta;
#if DUMMY != 0
    char dummy[DUMMY]; //32
#endif
    Branch branch[MAXCARD]; //32*55 = 1760 = 27.5  //1792 = 28
//#ifdef SHARED
    mutable std::shared_mutex* mutex_;
//#else
   // mutable std::shared_mutex* mutex_;
//#endif
};

#ifndef FULLLOG
struct splitLog{
    Node* parent; //+8
    Node* child; //+8
    Node* sibling; //+8
#ifdef MULTIMETA
    MetaData<META> meta;
#endif
};
#endif
#ifdef FULLLOG
struct splitLog{
    Node* parentPoint; //+8
    Node* childPoint; //+8
    Node* siblingPoint; //+8
    Node parent; 
    Node child; 
    Node sibling; 
#ifdef MULTIMETA
    MetaData<META> meta;
#endif
};
#endif

#ifndef INPLACE 
struct cowLog{
  MetaData<META> meta;
  Node* parent;
};
#endif


//#############################
#define METHODS 1
struct PartitionVars
{
    int partition[MAXCARD+1];
    int total, minfill;
    int taken[MAXCARD+1];
    int count[2];
    Rect cover[2];
    float area[2];
};

struct forSplit
{
    Branch BranchBuf[MAXCARD+1];
    int BranchCount;
    Rect CoverSplit;

    PartitionVars Partitions[METHODS];
};

/*
 * If passed to a tree search, this callback function will be called
 * with the ID of each data rect that overlaps the search rect
 * plus whatever user specific pointer was passed to the search.
 * It can terminate the search early by returning 0 in which case
 * the search will return the number of hits found up to that point.
 */

 //
 //extern struct Node * RTreeNewIndex();
 //extern struct Node * RTreeNewNode();
 //extern void RTreeInitNode(struct Node*);
 //extern struct Rect RTreeNodeCover(struct Node *);
 //extern void RTreeInitRect(struct Rect*);
 //extern struct Rect RTreeCombineRect(struct Rect*, struct Rect*);
 //extern int RTreeOverlap(struct Rect*, struct Rect*);
extern int RTreeAddBranch(PMEMobjpool *, struct Branch *, TOID(Node), TOID(Node)*, TOID(Node), TOID(splitLog));
 //extern int RTreePickBranch(struct Rect *, struct Node *);
 //extern void RTreeSplitNode(struct Node*, struct Branch*, struct Node*, struct Node*);
 #endif /* _INDEX_ */
