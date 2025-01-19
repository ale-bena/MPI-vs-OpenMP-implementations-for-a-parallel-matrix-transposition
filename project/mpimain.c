#include "mpifunc.h"
#define MAX 128

int main(int argc, char** argv){
//get arguments
	if(argc<3){
			printf("Insert a number!!");
			exit(1);
	}
  int n = atoi(argv[1]);
  if(isPowerOfTwo(n)==0){
      exit(1);
  }
  int tries=atoi(argv[2]);
  if(tries<1){
      exit(1);
  }

//init MPI environment
    int rank, size;
    srand(25);
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

// Ensure the number of processes divides n evenly
    if(size>n || n%size!=0){
      if(rank==0){
        printf("can't run because of size%d>n%d or n not multiple", size, n);
      }
      MPI_Finalize();
      return 0;
    }
//variables
    double s,e, time_elapsed, time_avg;
		double *sort;
    if(rank==0){
      sort = (double *)malloc(tries*sizeof(double));
		  if(sort==NULL){
			  printf("Memory allocation failed!!\n");
			  free(sort);
			  exit(1);
	    }
    }
    
//EXECUTION
 	float** M = NULL;//declareMatrix(n);
	float** T = NULL;//declareMatrix(n);
for(int i=0; i<tries; i++){
  cache_flush();
	//float** M = NULL;//declareMatrix(n);
	//float** T = NULL;//declareMatrix(n);
  M=declareMatrix(n);
  T=declareMatrix(n);
	if (rank == 0) { 
       initializeMatrix(M, n); 
       //printMatrix(M,n);
       s = MPI_Wtime();
    }
  matTransposeMPI(M, T, n, rank, size);
	if(rank==0){
		e = MPI_Wtime();
    //printMatrix(T,n);
		const char *times_table="times_table.csv";
		FILE *times = fopen(times_table, "a");
		if (times == NULL) {
			printf("Error opening file");
		exit(1); 
    }
		time_elapsed = (e-s);
		sort[i]=time_elapsed;
		fprintf(times, "MPI,%d,%d,%.9f\n",n, size, time_elapsed);
     fclose(times);
    }
  freeMatrix(M);
  freeMatrix(T);
  M=NULL; T=NULL;
}
//print on avg file
if(rank==0){
    bubbleSort(sort, tries);
    time_avg = timeAverage(sort, tries);
  	
     const char *times_average="times_average.csv";
	FILE *avg = fopen(times_average, "a+");
	if (avg == NULL) {
		printf("Error opening file");
		return 1;
  }
  
	char buffer[MAX];
	char string[] = "Sequential";
	double gettime = 0.0;
	fgets(buffer, MAX, avg);
	while(fgets(buffer, MAX, avg)){
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
	if(gettime==0.0){ //handles the case where the sequential program is not run before this program
     fprintf(avg, "MPI,%d,%d,%.9f,-,-\n", n, size, time_avg);
   }
   else{ 
    speed = speedup(gettime, time_avg); 
    eff = efficiency(speed, size);
	  fprintf(avg, "MPI,%d,%d,%.9f,%.9f,%.9f\n", n, size, time_avg, speed, eff);
	}
   fclose(avg);
   free(sort); 
   
 }
    MPI_Finalize();
    return 0;
}
