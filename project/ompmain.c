#include "ompfunc.h"
#define MAX 128

int main(int argc, char** argv){
	if(argc<2){
		printf("Insert a number!!");
		exit(1);
	}
	int n=atoi(argv[1]);
	if(isPowerOfTwo(n)==0){
		exit(1);
	}
	
	int n_threads=atoi(argv[2]);
	if(n_threads<1){
		exit(1);
	}
	omp_set_num_threads(n_threads);
	
	int tries=atoi(argv[3]);
	if(tries<1){
		exit(1);
	}
	
	srand(25);
	double s,e,time_elapsed;
	const char *times_table="times_table.csv";
	const char *header1 = "Type, Dim, Threads, Time";
	FILE *times = fopen(times_table, "a");
	if (times == NULL) {
		perror("Error opening file");
		return 1;
	}

//SORT ARRAY INITIALIZATION
	double *sort = (double *)malloc(tries*sizeof(double));
	if(sort==NULL){
		printf("Memory allocation failed!!\n");
		free(sort);
		exit(1);
	}
//EXECUTION
  float** M=NULL;
  float** T=NULL;
	for (int i=0; i<tries; i++){
		cache_flush();
		M=declareMatrix(n);
		initializeMatrix(M, n);
		T=declareMatrix(n);
		s=omp_get_wtime();
		matTransposeOMP(M,T,n);
		e=omp_get_wtime();
     if(n<=32){
      checkCorrectness(M,T,n);
    }
		time_elapsed = (e-s);
		sort[i]=time_elapsed;
		fprintf(times, "Omp,%d,%d,%.9f\n", n, n_threads, time_elapsed);		
		freeMemory(M);
    freeMemory(T);
    M=NULL; T=NULL;
    //freeMemory(M, n);
		//freeMemory(T, n);
	} 
	fclose(times);
 	bubbleSort(sort, tries);
 	double time_avg = timeAverage(sort, tries); 
  
 	const char *times_average="times_average.csv";
	FILE *avg = fopen(times_average, "a+"); 
	if (avg == NULL) {
		perror("Error opening file");
    exit(1);
	}

	char buffer[MAX];
	char string[] = "Sequential";
	double gettime = 0.0;
	fgets(buffer, MAX, avg);
	while(fgets(buffer, MAX, avg)!=NULL){
		char *token = strtok(buffer, ",");
		int i=0;
		int flag=0;
		while(token != NULL){
			i++;
			if(i==1){
				if(strcmp(token, string)==0){
					flag=1;
					//printf("\nflag 1.1\n");
				}
			}
			if(i==2){
				int tmp = atoi(token);
				if(tmp==n && flag==1){
					flag=1;
					//printf("flag 1.2\n");
				}
				else{
					flag=0;
				} 
			}
			
			if(i==4 && flag==1){
				gettime=atof(token);
			}
			token = strtok(NULL, ",");
		}
	
	}
	double speed=0.0;
	double eff=0.0;
	//printf("\nGettime:%.9f\nAvgTime:%.9f\nSpeedup:%.9f\nEfficiency:%.9f", gettime, time_avg, speed, eff);
   if(gettime==0.0){  //handles the case where the sequential program is not run before this program
     fprintf(avg, "Omp,%d,%d,%.9f,-,-\n", n, n_threads, time_avg);
   }
   else{
	  speed = speedup(gettime, time_avg);
	  eff = efficiency(speed, n_threads);
	  fprintf(avg, "Omp,%d,%d,%.9f,%.9f,%.9f\n", n, n_threads, time_avg, speed, eff);
   }
	fclose(avg);
  free(sort); 
  return 0;
}
