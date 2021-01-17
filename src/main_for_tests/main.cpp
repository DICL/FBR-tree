#ifdef INPLACE
#include "../FBRtree_INP.h"
#else
#include "../FBRtree_COW.h"
#endif
#include <unistd.h>
#include <time.h>
#include <stdlib.h>

#ifdef BREAKDOWN
extern double traversal_time;
extern double write_time;
double k_insert_time = 0.0;
double k_insert_time_v2 = 0.0;
#ifdef SPLIT
extern double split_write_time;
double k_split_insert_time = 0.0;
double k_split_insert_time_v2 = 0.0;
#endif
#endif
extern int isSplit;
uint64_t clflushCnt = 0;
uint64_t flipCount = 0;
uint64_t k_tmp_clflushCnt = 0;
uint64_t k_tmp_flipCount = 0;
#ifdef SPLIT
extern int split_cnt;
uint64_t k_split_clflushCnt = 0;
uint64_t k_split_flipCount = 0;
#endif

double conc_time = 0.0;
double search_time = 0.0;
double insert_time = 0.0;

#define POOL_SIZE (10737418240)

char queryFileName[] = "/home/bnam/FBRtree/input1M.txt";
char dataFileName[] = "/home/bnam/FBRtree/input1M.txt";
int NUMDATA=80000000;
int SEARCH=10000;
int IpthreadNum=8;

void cache_clear(){
  // Remove cache
  uint64_t size = 1024*1024*1024;
  char *garbage = new char[size];
  for (int i=0; i<size; ++i){
      garbage [i] = i;
  }
  for (int i=100; i<size; ++i){
      garbage [i] += garbage[i-100];
  }
  
  delete[] garbage;
}

TOID(Node) total_root = OID_NULL;

// thread data
struct thread_data{
  PMEMobjpool *pop;
  int  tid =0;
  uint64_t  hit=0;
  struct Rect* rect=NULL;
  int startDataNum=0; // start number of fraction for each thread
  int dataNum=0;      // number of data for each thread
};

#ifdef CONC
struct thread_data_both{
  	PMEMobjpool *pop;
	int tid=0;
	int hit=0;
	Rect* Irect = NULL;
	Rect* Srect = NULL;
	int startDataNum=0;	//start number of fraction for each thread
	int dataINum=0;		//number of data for each thread
	int dataSNum=0;		//number of data for each thread
};

void* PThreadBOTH(void *arg){
	Rect r;
	thread_data_both* td = (struct thread_data_both*) arg;

 	struct timespec begin,end;
	int iCur=0;
	int sCur=0;
	for(int i=0; i< td->dataINum + td->dataSNum; i++){
		clock_gettime(CLOCK_REALTIME, &begin);
		if(rand()%10 <1){
			RTreeInsertRect(td->pop, &(td->Irect[iCur]), td->startDataNum+iCur+1, total_root, OID_NULL);
			iCur++;
		}else{
			td->hit += hostRTreeSearch(total_root, &(td->Srect[sCur]));
			sCur++;
		}
		clock_gettime(CLOCK_REALTIME, &end);
	
		double dur = (end.tv_sec - begin.tv_sec)*1000*1000*1000 + (end.tv_nsec - begin.tv_nsec); 

    conc_time += dur;
 }

 pthread_exit(NULL);
}
#endif

void* PThreadInsert(void *arg)
{
  // Insert the R-tree boundary region
  struct Rect r;
  struct thread_data* td = (struct thread_data*) arg;
  struct timespec begin,end;

  for(int i=0; i<td->dataNum; i++){
	  isSplit = 0;
	  k_tmp_clflushCnt = 0;
	  k_tmp_flipCount = 0;
	  
    clock_gettime(CLOCK_REALTIME, &begin);
    RTreeInsertRect(td->pop, &(td->rect[i]), td->startDataNum+i+1, (total_root), OID_NULL);
	  clock_gettime(CLOCK_REALTIME, &end);

	  double dur = (end.tv_sec - begin.tv_sec)*1000*1000*1000 + (end.tv_nsec - begin.tv_nsec); 
	  insert_time += dur;

	  clflushCnt += k_tmp_clflushCnt;
	  flipCount += k_tmp_flipCount;
#ifdef BREAKDOWN
    k_insert_time_v2 += dur;
#endif 

#ifdef SPLIT
    if (isSplit) {
      split_cnt++;
#ifdef BREAKDOWN
      k_split_insert_time_v2 += dur;
#endif
      k_split_clflushCnt += k_tmp_clflushCnt;
      k_split_flipCount += k_tmp_flipCount;
    }
#endif
  }
  pthread_exit(NULL);
}


// P Thread Search
void* PThreadSearch(void* arg)
{
  int i;
  struct Rect r;
  struct thread_data* td = (struct thread_data*) arg;

  struct timespec begin,end;
  
  
  for(i=0;i<td->dataNum;i++){
  	clock_gettime(CLOCK_REALTIME, &begin);
    td->hit += hostRTreeSearch(total_root, &(td->rect[i]));
	  clock_gettime(CLOCK_REALTIME, &end);
  
    double dur = (end.tv_sec - begin.tv_sec)*1000*1000*1000 + (end.tv_nsec - begin.tv_nsec); 
 	  search_time += dur; 
  }
  pthread_exit(NULL);
}

Rect *r = NULL;
Rect *sr = NULL;

//-----------------------------------------------------------------------
//------------------------------- MAIN ----------------------------------
//-----------------------------------------------------------------------
int main(int argc, char *args[])
{
	// Check the arguments for extra information
	if(argc<6){
#ifdef CONC
        printf("Usage: %s path (number_of_INSERT) (number_of_SEARCH) (number_of_TOTAL_THREADs) (write_Latency) (k_delta)\n", args[0]);
#else
        printf("Usage: %s path (number_of_INSERT) (number_of_SEARCH) (number_of_Insert_THREADs) (number_of_Search_THREADs) (write_Latency) (k_delta)\n", args[0]);
#endif
	    exit(1);
	}
	
    NUMDATA = atoi(args[2]);	// Initialize the number of Data
    SEARCH = atoi(args[3]);	// Initialize the number of search Data
    IpthreadNum = atoi(args[4]);	// Initialize the number of insert Thread
#ifdef CONC
    writeLatency = 0;	// Initialize the number of Thread 
    float delta = atof(args[6]);
    printf("INSERT: %d, SEARCH: %d, thread: %d, Write_Latency: %d delta: %f * FLT_MIN\n", NUMDATA, SEARCH, IpthreadNum, writeLatency, delta);
  	delta = delta * FLT_MIN;
#else
	  int SpthreadNum = atoi(args[5]);	// Initialize the number of search Thread     
    writeLatency = atoi(args[6]);	// Initialize the number of Thread 
    float delta = atof(args[7]);
    printf("INSERT: %d, SEARCH: %d, insert_thread: %d, search_thread: %d, Write_Latency: %d delta: %f * FLT_MIN\n", NUMDATA, SEARCH, IpthreadNum, SpthreadNum, writeLatency, delta);
#endif
    printf("=== NODECARD=%d\n",NODECARD);

	//##################### GET DATA #########################
    FILE* fp = fopen(dataFileName, "r+b");

    if(fp==0x0){
        printf("Line %d : Insert file open error\n", __LINE__);
        exit(1);
    }

    r = new Rect[NUMDATA];
    for(int i=0; i<NUMDATA; i++){
        fscanf(fp, "%f %f %f %f %f %f", &r[i].boundary[0], &r[i].boundary[1], &r[i].boundary[2], &r[i].boundary[3], &r[i].boundary[4], &r[i].boundary[5]);
   }
	
    if(fclose(fp) != 0) printf("Insert file close error\n");

    //################ GET SEARCH DATA #######################
    fp = fopen(queryFileName, "r+b");

    if(fp==0x0){
        printf("Line %d : Search file open error\n", __LINE__);
        exit(1);
    }

    sr = new Rect[SEARCH];
    for(int i=0; i<SEARCH; i++){
      fscanf(fp, "%f %f %f %f %f %f", &sr[i].boundary[0], &sr[i].boundary[1], &sr[i].boundary[2], &sr[i].boundary[3], &sr[i].boundary[4], &sr[i].boundary[5]);

      sr[i].boundary[0] -= delta;
      sr[i].boundary[1] -= delta;
      sr[i].boundary[2] -= delta;
      sr[i].boundary[3] += delta;
      sr[i].boundary[4] += delta;
      sr[i].boundary[5] += delta;
    }
    if(fclose(fp) != 0) 
        printf("Search file close error\n");
	
    printf("FInish prepare DATA\n");
    //########################################################

    //################ PREPARE TO PM ######################
    char path[32];
    strcpy(path, args[1]);

    PMEMobjpool *pop;
    if(access(path, 0) != 0){ 
        pop = pmemobj_create(path, "rtree", POOL_SIZE, 0666);
        if(pop == NULL){
    			perror("pmemobj_create");
		    return 1;
		}
		total_root = POBJ_ROOT(pop, Node);
		char *ptr = (char*)((unsigned long)(char*)total_root.oid.off &~(cacheLineSize-1)); 
		total_root.oid.off = (uint64_t)ptr; 

		RTreeNewIndex(pop, total_root);	
	}else{
		pop = pmemobj_open(path, "rtree");
		total_root = POBJ_ROOT(pop, Node);
		char *ptr = (char*)((unsigned long)(char*)total_root.oid.off &~(cacheLineSize-1)); 
		total_root.oid.off = (uint64_t)ptr; 
		if(total_root.oid.off){    
			RTreePrint(total_root);
			return 1;
		}
  }
  printf("FInish prepare PM\n");
	//########################################################
	//################ CREATE & WARM_UP ######################
  // CREATE RTree total_root node
#ifdef CONC
	// insert : search = 1 : 9
	const int num_insert_mixed = 2000;
	const int num_search_mixed = 18000;
	const int warmData = NUMDATA - num_insert_mixed;
	NUMDATA = num_insert_mixed;
	SEARCH = num_search_mixed;
  printf("MIXED_INSERT: %d, MIXED_SEARCH: %d\n", NUMDATA, SEARCH);
	
	for(int i=0; i<warmData; i++){
		RTreeInsertRect(pop, &r[i], i+1, total_root, OID_NULL);
	}
	printf("[WARM UP] end %d/%d\n", warmData, NUMDATA);
#else
    const int insertPerThread = NUMDATA / IpthreadNum;
    const int searchPerThread = SEARCH / (SpthreadNum);
#endif

#ifndef CONC
	//################# MULTI INSERT #########################
    struct timeval it1,it2;
    double time_it1;

    printf("[Only Concurrent inserting]\n");
    int rc;
    void *status;
    pthread_t *threads_i = new pthread_t[IpthreadNum];
    struct thread_data* td_i = new thread_data[IpthreadNum];
    
    gettimeofday(&it1,0); // start the stopwatch
    for(int i=0; i<IpthreadNum; i++){
        int dataS = i*insertPerThread;
        int dataE = (i == IpthreadNum-1)? NUMDATA-dataS : insertPerThread;
       
		    td_i[i].pop = pop;
        td_i[i].tid = i;
        td_i[i].hit = 0;
        td_i[i].rect = &r[dataS];
        td_i[i].startDataNum = dataS;
        td_i[i].dataNum = dataE;
    
        rc = pthread_create(&threads_i[i], NULL, PThreadInsert, (void *)&td_i[i]);
        if (rc){
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1); 
		    }
	  }

    for(int i=0; i<IpthreadNum; i++){
        rc = pthread_join(threads_i[i], &status);
        if (rc) {
            printf("ERROR; return code from pthread_join() is %d\n", rc);
            exit(-1);
        }
    }
  	gettimeofday(&it2,0);
    time_it1 = (it2.tv_sec-it1.tv_sec)*1000000 + (it2.tv_usec - it1.tv_usec);
   
#ifdef BREAKDOWN
	printf("throughput        : %.3lf\n", NUMDATA/(time_it1));
	printf("TOTAL insert time (usec)        : %lf\n", k_insert_time_v2);
	printf("AVG insert time (usec)          : %.3lf\n", k_insert_time_v2 / NUMDATA);
	printf("AVG TRAVERSAL TIME (usec/insertion)              : %.3lf\n", traversal_time/NUMDATA);

#ifdef SPLIT
	printf("num_split_insertion: %d\n", split_cnt);
	printf("TOTAL split insert time (usec)  : %lf\n", k_split_insert_time_v2);
	printf("AVG split insert time (usec/insertion): %lf\n", k_split_insert_time_v2/split_cnt);
#endif

  printf("clflush: %ld, flipCount: %ld\n", clflushCnt, flipCount); 	
	printf("AVG clflush (per insertion): %.3lf\n", (double) clflushCnt / NUMDATA);
	printf("AVG flipCount (per insertion): %.3lf\n", (double) flipCount / NUMDATA);
#ifdef SPLIT
  printf("TOTAL split_clflush: %ld, TOTAL split_bitflip: %ld\n", k_split_clflushCnt, k_split_flipCount);
	printf("AVG SPLIT CLFLUSH COUNT (per insertion)      : %.3lf\n", (double) k_split_clflushCnt / split_cnt);
	printf("AVG SPLIT BITFLIP COUNT (per insertion)          : %.3lf\n", (double) k_split_flipCount / split_cnt);
#endif

#else
	printf("throughput        : %lf\n", NUMDATA/(insert_time));
	printf("TOTAL insert time (nsec)        : %lf\n", insert_time);
	printf("AVG insert time (usec)          : %lf\n", insert_time / NUMDATA);
	printf("\n");
	double total_insert_time_milli = time_it1 / 1000;
	printf("TOTAL insert time + extra (msec)    : %lf\n", total_insert_time_milli);
	printf("AVG insert time + extra (msec)      : %lf\n", total_insert_time_milli / NUMDATA);
	printf("INSERT THROUGHPUT (insert/msec) : %lf\n", NUMDATA / total_insert_time_milli);
#endif

  printf("\n");
	fprintf(stderr, "Insertion is done.\n");
  cache_clear();

	//########################################################
//-------------------------------------------------------------------------
   
	//###################  SEARCH DATA #######################
  struct timeval t1,t2;
  double time_t2;
  uint64_t hit = 0;

	printf("[Only Concurrent Searching] %d\n", SpthreadNum);
	pthread_t* threads = new pthread_t[SpthreadNum];
	struct thread_data* td = new thread_data[SpthreadNum];

	gettimeofday(&t1,0); // start the stopwatch
	for(int i=0; i<SpthreadNum; i++){
		int searchS = i*searchPerThread;
		int searchE = (i==SpthreadNum-1)? SEARCH-searchS : searchPerThread;

		td[i].tid = i;
		td[i].hit = 0;
		td[i].rect = &sr[searchS];
		td[i].dataNum = searchE;

		rc = pthread_create(&threads[i], NULL, PThreadSearch, (void *)&td[i]);
		if (rc){
			printf("ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}
	}

	for(int i=0; i<SpthreadNum; i++){
		rc = pthread_join(threads[i], &status);
		if (rc) {
			printf("ERROR; return code from pthread_join() is %d\n", rc);
			exit(-1);
		}
		hit += td[i].hit;
	}
	gettimeofday(&t2,0);
	time_t2 = (t2.tv_sec-t1.tv_sec)*1000000 + (t2.tv_usec - t1.tv_usec);

	
	printf("host search time (nsec): %lf\n", search_time);
	printf("throughput             : %lf\n", SEARCH/(search_time));
	printf("Host Hit counter = %ld, PThread: %d\n", hit, SpthreadNum);

	printf("single thread: AVG search time(usec): %lf\n", search_time/SEARCH);

	printf("\n");
	double total_search_time_milli = time_t2 / 1000;
	printf("TOTAL search time + extra (msec)    : %lf\n", total_search_time_milli);
	printf("AVG search time + extra (msec)      : %lf\n", total_search_time_milli / SEARCH);
	printf("SEARCH THROUGHPUT (search/msec) : %lf\n", SEARCH / total_search_time_milli);
	
  hit = 0;

  printf("done-------------------------------------------------------\n");

#else
	printf("[CONCURRENT!!!]\n");
	struct timeval t1, t2;
	double time_t2;
	int rc;
	void* status;
	uint64_t hit = 0;

	int insertThreads = IpthreadNum;
	int searchThreads = IpthreadNum;

	pthread_t threads[IpthreadNum];
	struct thread_data_both td[IpthreadNum];

	int in=0, se=0;

	const int insertPerThread = (NUMDATA) / insertThreads;
	const int searchPerThread = SEARCH / searchThreads;

	gettimeofday(&t1, 0);
	for(int i=0; i<IpthreadNum; i++){
		int dataS = warmData + in*insertPerThread;
		int dataE = (in == insertThreads-1)? NUMDATA-in*insertPerThread : insertPerThread;
		in++;
		int searchS = se*searchPerThread;
		int searchE = (se== searchThreads-1)? SEARCH-searchS : searchPerThread; 
		se++;

		td[i].pop = pop;
		td[i].tid = i;
		td[i].hit = 0;
		td[i].Irect = &r[dataS]; 
		td[i].Srect = &sr[searchS];
		td[i].startDataNum = dataS;
		td[i].dataINum = dataE;
		td[i].dataSNum = searchE;

		rc = pthread_create(&threads[i], NULL, PThreadBOTH, (void *)&td[i]);
		if (rc) {
			printf("ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}
	}

	for(int i=0; i<IpthreadNum; i++){
		rc = pthread_join(threads[i], &status);
		if (rc) {
			printf("ERROR; return code from pthread_join() is %d\n", rc);
			exit(-1);
		}
		hit += td[i].hit;
	}

	gettimeofday(&t2,0);
	time_t2 = (t2.tv_sec-t1.tv_sec)*1000000 + (t2.tv_usec - t1.tv_usec);

	printf("TOTAL concurrent time (nsec)        : %lf\n", conc_time);
	printf("throughput            : %lf\n", (20000)/conc_time);
	printf("Host Hit counter = %ld\n", hit);
	printf("\n");
	
	printf("AVG concurrent time (nsec)  : %lf\n", conc_time/ (SEARCH+NUMDATA));
	printf("AVG search  time (nsec)          : %lf\n", conc_time / SEARCH);
	printf("\n");

	double total_time_milli = time_t2 / 1000;
	printf("TOTAL execution time + extra (msec) : %lf\n", total_time_milli);
	printf("MIXED THROUGHPUT (query/msec)   : %lf\n", (NUMDATA + SEARCH) / total_time_milli);
	printf("\n");
#endif

	pmemobj_close(pop);
	return 0;
}
