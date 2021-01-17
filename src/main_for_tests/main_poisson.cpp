// *************************
// This is test for poisson.
// We insert 200 data during search 1048576 data with 48 threads
// ./poissonTest {pmem_path} 200 1048576 1 48 0 16082
// *************************
#include "../FBRtree_poi.h"
#include "../poisson.h"
#include <unordered_map>
#include <time.h>
#include <unistd.h>
#include <cfloat>
#include <solar.h>

int NUMDATA=100;
int SEARCH=10;
int IpthreadNum=8;

//poisson
class Query{
	public:
		Rect qr;
		unsigned long long arrival_time;
};
extern unsigned long long poisson(double x);
//

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

//generate search data
Rect* sr = NULL;
int selection_ratio = 1;
void makeSearch(int num){
	float f_r = selection_ratio;
	sr = new Rect[num];
	srand(0);
	for(int i=0; i<num; i++){
		sr[i].boundary[0] = 0.0;
		sr[i].boundary[2] = 0.0;
		sr[i].boundary[3] = FLT_MAX;
		sr[i].boundary[5] = FLT_MAX;
		float a = (float)(rand() % (100 - (selection_ratio-1)));
		sr[i].boundary[1] = a + 0.49;
		sr[i].boundary[4] = a + (0.01+f_r);
	}
}


TOID(Node) total_root = OID_NULL;

// thread data
struct thread_data{
  PMEMobjpool *pop;
  int  tid =0;
  uint64_t  hit=0;
  struct kv_t* InsKV=NULL;
  char* filePath = NULL;
  int startDataNum=0; // start number of fraction for each thread
  int dataNum=0;      // number of data for each thread
  queue<Query*>* que = NULL;
  vector<unsigned long long>* lv = NULL;
};

bool ended = false;
void* PThreadInsert(void *arg)
{
  // Insert the R-tree boundary region
  sleep(180);
  struct thread_data* td = (struct thread_data*) arg;
  int j = 0;
  for(int i=0; i<100; i++){
	TOID(splitLog) slog = total_log;
	slog.oid.off = total_log.oid.off + sizeof(splitLog)*(td->tid*2);
    RTreeInsertRect(td->pop, (Rect*)&(td->InsKV[i].key), td->filePath, (total_root), slog);
  }

  pthread_exit(NULL);
}

// P Thread Search
Rect *r = NULL;

unsigned long long avg_time= 0;
void* PThreadSearch(void* arg)
{
	struct thread_data* td = (struct thread_data*) arg;
	struct timespec s, cur, f;
	unsigned long long cur_time= 0;
	clockid_t clk_id = CLOCK_MONOTONIC;
	clock_gettime(clk_id, &s);
	clock_gettime(clk_id, &cur);
	cur_time = (cur.tv_sec - s.tv_sec)*1000000000LLU + (cur.tv_nsec - s.tv_nsec);
	for(int i=0; i < td->dataNum; i++){
		int reSplit = 1;
		Query* q = td->que->front();
		td->que->pop();
		while(q->arrival_time > cur_time){
			clock_gettime(clk_id, &cur);
			cur_time = (cur.tv_sec - s.tv_sec)*1000000000LLU + (cur.tv_nsec - s.tv_nsec);
		}
		
		td->hit += hostRTreeSearch(total_root, &(q->qr));

		clock_gettime(clk_id, &cur);
		cur_time = (cur.tv_sec - s.tv_sec)*1000000000LLU + (cur.tv_nsec - s.tv_nsec);
		td->lv->push_back(cur_time - q->arrival_time);
	}
	clock_gettime(clk_id, &f);

  pthread_exit(NULL);
}


//-----------------------------------------------------------------------
//------------------------------- MAIN ----------------------------------
//-----------------------------------------------------------------------
int main(int argc, char *args[])
{
	// Check the arguments for extra information
	if(argc<6){
        printf("Usage: %s path (number_of_INSERT) (number_of_SEARCH) (number_of_Insert_THREADs) (number_of_Search_THREADs) (write_Latency) (fileN) \n", args[0]);
   	    exit(1);
	}
	
    NUMDATA = atoi(args[2]);	// Initialize the number of Data
    SEARCH = atoi(args[3]);	// Initialize the number of search Data
    IpthreadNum = atoi(args[4]);	// Initialize the number of insert Thread
	  int SpthreadNum = atoi(args[5]);	// Initialize the number of search Thread     
    writeLatency = atoi(args[6]);	// Initialize the number of Thread 
    int fileN = atoi(args[7]);
    selection_ratio = 1.0; 
	  printf("INSERT: %d, SEARCH: %d, insert_thread: %d, search_thread: %d, Write_Latency: %d fileN: %d, selection_ratio: %d\n", NUMDATA, SEARCH, IpthreadNum, SpthreadNum, writeLatency, fileN, selection_ratio);
    
	printf("NODECARD=%d\n",NODECARD);

	//##################### GET DATA #########################
	struct timeval f1,f2;
	double time_f2;
	int arr_size = 100 * fileN; //warmup;                         
	kv_t* kv_arr = new kv_t[arr_size];  

	gettimeofday(&f1,0); 
	get_kv_narr(fileN, kv_arr); 
	gettimeofday(&f2,0);
	time_f2 = (f2.tv_sec-f1.tv_sec)*1000000 + (f2.tv_usec - f1.tv_usec);  
	printf("HDF5 file Num        : %d\n", fileN); 
	printf("HDF5 file open time  : %f\n", time_f2); 
	printf("throughput           : %.3lf\n", (1000)/(time_f2));  
  //################ GET SEARCH DATA #######################
 	makeSearch(SEARCH);	
  //########################################################

  //################ PREPARE TO PM ######################
    char path[32];
    strcpy(path, args[1]);

    PMEMobjpool *pop;
    if(access(path, 0) != 0){ 
        pop = pmemobj_create(path, "rtree", 1024*1024*1024, 0666);
        if(pop == NULL){
			perror("pmemobj_create");
		    return 1;
		}
		total_root = POBJ_ROOT(pop, Node);
		RTreeNewIndex(pop, total_root);	
		log_init(pop, IpthreadNum);
   
	}else{
		pop = pmemobj_open(path, "rtree"); 
		total_root = POBJ_ROOT(pop, Node);
		if(total_root.oid.off){    
			RTreePrint(total_root);
			return 1;
		}
  }
	//########################################################
	
  //################ CREATE & WARM_UP ######################
	const int warmData = (fileN-2)*100; // warm-up almost data except only 200 data
	for(int i=0; i< warmData; i++){
		RTreeInsertRect(pop, (struct Rect*)&kv_arr[i].key, kv_arr[i].val, total_root, OID_NULL); 
	}
	printf("Warm up (warm of NUMDATA) end %d/%d\n", warmData, NUMDATA );

	const int insertPerThread = (NUMDATA) / IpthreadNum;  // 200 / 2
	const int searchPerThread = SEARCH / SpthreadNum; // 48*1M / 48 = 1M
	printf("insert: %d, search: %d\n", insertPerThread, searchPerThread);

	//################# INSERT #########################
  int rc;
  void *status;
  pthread_t *threads_i = new pthread_t[IpthreadNum];
  struct thread_data* td_i = new thread_data[IpthreadNum];
   
	for(int i=0; i<IpthreadNum; i++){
    int dataS = warmData;    
    int dataE = 100;

		td_i[i].pop = pop;
		td_i[i].tid = i;
		td_i[i].hit = 0;
		td_i[i].InsKV = &kv_arr[dataS];  
		td_i[i].filePath = kv_arr[dataS].val; 
		td_i[i].startDataNum = dataS;
		td_i[i].dataNum = dataE;
	}
	fprintf(stderr, "Insertion is done.\n");

	cache_clear();

	//########################################################
	queue<Query*>* qq = new queue<Query*>[SpthreadNum];
	unsigned long long arrival_time=0;
	rand_val(1);
	srand(1);
	vector< unsigned long long>* latency_vector = new vector< unsigned long long>[SpthreadNum];

	//###################  SEARCH DATA #######################
	struct timeval t1,t2;
	double time_t2;
	uint64_t hit = 0;

	printf("[Concurrent Searching]\n");

	pthread_t threads[SpthreadNum];
	struct thread_data td[SpthreadNum];

	int j=0,  h=0; 
	gettimeofday(&t1,0); // start the stopwatch
	for(int i=0; i<SpthreadNum; i++){
		int searchS = i*searchPerThread;
		int searchE = (i==SpthreadNum-1)? SEARCH-searchS : searchPerThread;
		unsigned long long old_arrival = 0;
		arrival_time = 0;
		for(int g=0; g< searchE; g++){

#ifdef SHARED
			arrival_time += poisson(1.0/40786.098)*10000;
#else
				arrival_time += poisson(1.0/26972.003)*10000;
#endif
			
      Query *q = new Query;
			q->qr = sr[searchS + g];
			q->arrival_time = arrival_time ;
			qq[i].push(q);
		}
		td[i].que = &qq[i];
		td[i].lv = &latency_vector[i];

		td[i].pop = pop;
		td[i].tid = i;
		td[i].hit = 0;
		td[i].startDataNum = searchS;
		td[i].dataNum = searchE;

		rc = pthread_create(&threads[i], NULL, PThreadSearch, (void *)&td[i]);
		if (rc){
			printf("ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}
	}
	rc = pthread_create(&threads_i[0], NULL, PThreadInsert, (void *)&td_i[0]);

	rc = pthread_join(threads_i[0], &status);
	if (rc) {
		printf("ERROR; return code from pthread_join() is %d\n", rc);
		exit(-1);
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
	printf("Host total time (msec): %.3lf\n", time_t2/1000);
	printf("Host Hit counter = %ld, PThread: %d\n", hit, SpthreadNum);
	//########################################################

	pmemobj_close(pop);
	return 0;
}
