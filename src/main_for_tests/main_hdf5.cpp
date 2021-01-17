#include "../FBRtree_poi.h"
#include <unistd.h>
#include <cfloat>
#include <solar.h>

extern uint64_t flipCount;
extern uint64_t splitCount;
int NUMDATA=100;
int SEARCH=10;
int IpthreadNum=8;

void cache_clear() {
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

// thread data
struct thread_data {
  PMEMobjpool *pop; 
  int  tid =0;
  uint64_t  hit=0;
  struct Rect* rect=NULL;
  char* filePath = NULL;
  int startDataNum=0; // start number of fraction for each thread
  int dataNum=0;      // number of data for each thread
};

TOID(Node) total_root = OID_NULL;
//P Thread Insert
#ifdef CONC
struct thread_data_both {
  PMEMobjpool *pop;
  int  tid =0;
  int  hit=0;
  struct kv_t* InsKV=NULL;
  struct Rect* SerKV=NULL;
  char* filePath = NULL;
  int startDataNum=0; // start number of fraction for each thread
  int dataINum=0;      // number of data for each thread
  int dataSNum=0;      // number of data for each thread
};

struct Rect* srch = NULL;
void makeSearch(int num) {
  srch = new Rect[num];
  srand(0);
  for(int i=0; i< num; i++) {
      srch[i].boundary[0] = 0.0;
      srch[i].boundary[3] = FLT_MAX;
      srch[i].boundary[2] = 0.0;
      srch[i].boundary[5] = FLT_MAX;
      float a = (float)(rand() % 90);
      srch[i].boundary[1] = a;
      srch[i].boundary[4] = a + 9.49;
  }
}

void* PThreadBOTH(void *arg)
{
  // Insert the R-tree boundary region
  struct thread_data_both* td = (struct thread_data_both*) arg;

  int iCur = 0;
  int sCur = 0;
#ifdef OS
  for(int i=0; sCur < td->dataSNum; i++) {
    	  td->hit += hostRTreeSearch(total_root, &(td->SerKV[sCur]));
      	sCur++;
  }
#else
  for(int i=0; i < td->dataINum + td ->dataSNum; i++) {
	  if(iCur < td->dataINum){
		  TOID(splitLog) sLog = total_log;
		  sLog.oid.off = total_log.oid.off + sizeof(splitLog)*(td->tid*2);
          RTreeInsertRect(td->pop, (Rect*)&(td->InsKV[iCur].key), td->InsKV[iCur].val, total_root, sLog);
		  iCur++;
	  } else {
	      //TODO: 60% is Search
    	  td->hit += hostRTreeSearch(total_root, &(td->SerKV[sCur]));
      	sCur++;
	  }
  }
#endif
  pthread_exit(NULL);
}
#endif


//-----------------------------------------------------------------------
//------------------------------- MAIN ----------------------------------
//-----------------------------------------------------------------------
int main(int argc, char *args[])
{
  // Check the arguments for extra information
  if(argc<4){
    printf("Usage: %s path (number_of_INSERT) (number_of_SEARCH) (number_of_total_THREADs)\n", args[0]);
    exit(1);
  }

  NUMDATA = atoi(args[2]);	// Initialize the number of Data
  SEARCH = atoi(args[3]);	// Initialize the number of search Data
  IpthreadNum = atoi(args[4]);	// Initialize the number of insert Thread
  int fileN = NUMDATA/100;	// Initialize the write latency 
  printf("INSERT: %d, SEARCH: %d, total_thread: %d, fileN: %d \n", NUMDATA, SEARCH, IpthreadNum, fileN);

  struct timeval f1,f2;
  double time_f2;
  //##################### GET DATA #########################
#ifndef OS
  fileN = 16081; //100 data
#endif
  int arr_size = 100 * fileN;
  kv_t* kv_arr = new kv_t[arr_size];

  gettimeofday(&f1,0); // start the stopwatch
  get_kv_narr(fileN, kv_arr);
  gettimeofday(&f2,0); // start the stopwatch
  time_f2 = (f2.tv_sec-f1.tv_sec)*1000000 + (f2.tv_usec - f1.tv_usec);
  printf("HDF5 file Num        : %d\n", fileN);
  printf("Read HDF5 time (msec): %.3lf\n", time_f2/1000);
  printf("throughput           : %.3lf\n", (1000)/(time_f2));

  makeSearch(SEARCH);
  //################ PREPARE FOR PM  ######################
	char path[32];
	strcpy(path, args[1]);

	PMEMobjpool *pop;
	if(access(path, 0) != 0) {  
		pop = pmemobj_create(path, "rtree", 10737418240, 0666);
    if(pop == NULL) {
      perror("pmemobj_create");
			return 1;
		}
		total_root = POBJ_ROOT(pop, Node);
		char *ptr = (char*)((unsigned long)(char*)total_root.oid.off &~(cacheLineSize-1));
		total_root.oid.off = (uint64_t)ptr;

		RTreeNewIndex(pop, total_root);	
		log_init(pop, IpthreadNum);
	} else {
		pop = pmemobj_open(path, "rtree");
		total_root = POBJ_ROOT(pop, Node);
		char *ptr = (char*)((unsigned long)(char*)total_root.oid.off &~(cacheLineSize-1));
		total_root.oid.off = (uint64_t)ptr;
		if(total_root.oid.off) {    
			RTreePrint(total_root);
			return 1;
		}
	}

	//################ CREATE & WARM_UP ######################
    // CREATE RTree total_root node
    uint64_t ihit=0;

#ifdef OS
    const int warmData = NUMDATA; 
#else
    const int warmData = 500000; 
#endif

    for(int i=0; i< warmData; i++) {
        RTreeInsertRect(pop, (struct Rect*)&kv_arr[i].key, kv_arr[i].val, total_root, OID_NULL);
    }
    printf("Warm up (warm of NUMDATA) end %d/%d\n", warmData, NUMDATA);
    
    //########################################################
    struct timeval t1,t2;
    double time_t2;
    int rc;
    void *status;
    
    uint64_t hit = 0;

    int pthreads = IpthreadNum;

    pthread_t threads[pthreads];
    struct thread_data_both td[pthreads];
    
    int in=0, se=0; 

    const int insertPerThread = (NUMDATA) / pthreads;
    const int searchPerThread = SEARCH / pthreads;

    printf("Warm up %d (warm of NUMDATA) end %d/%d\n", warmData, NUMDATA);
    gettimeofday(&t1,0); // start the stopwatch
    
    for(int i=0; i<pthreads; i++) {
      int searchS = se*searchPerThread;
      int searchE = (se== pthreads-1)? SEARCH-searchS : searchPerThread; 
     
      se++;
	    td[i].pop = pop;
      td[i].tid = i;
      td[i].hit = 0;
      td[i].SerKV = &srch[searchS];
      td[i].dataSNum = searchE;

#ifndef OS      
      int dataS = warmData + in*insertPerThread;
      int dataE = (in == pthreads-1)? NUMDATA-in*insertPerThread : insertPerThread;
      in++;
      
      td[i].InsKV = &kv_arr[dataS];
      td[i].filePath = kv_arr[dataS].val;
      td[i].startDataNum = dataS;
      td[i].dataINum = dataE;
#endif
     
     rc = pthread_create(&threads[i], NULL, PThreadBOTH, (void *)&td[i]);
     if (rc) {
        printf("ERROR; return code from pthread_create() is %d\n", rc);
       exit(-1);
     }
   }

   for(int i=0; i<IpthreadNum; i++) {
      rc = pthread_join(threads[i], &status);
      if (rc) {
        printf("ERROR; return code from pthread_join() is %d\n", rc);
        exit(-1);
      }
      hit += td[i].hit;
   }

   gettimeofday(&t2,0);
   time_t2 = (t2.tv_sec-t1.tv_sec)*1000000 + (t2.tv_usec - t1.tv_usec);
   printf("concurrent time (msec): %.3lf\n", time_t2/1000);
   printf("Num of Insert & Search: %d\n", NUMDATA+SEARCH);
   printf("throughtput		  : %.3lf\n", (SEARCH)/(time_t2/1000));
   printf("Host Hit counter = %ld\n", hit);

   pmemobj_close(pop);
   return 0;
}
