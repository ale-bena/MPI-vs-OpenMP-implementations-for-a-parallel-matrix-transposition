#include "mpifunc.h"
#define MAX 128

int main(int argc, char** argv){
//get arguments: dim e thread(gia presi con MPI)
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
/*		const char *times_table="times_table.csv";
		FILE *times = fopen(times_table, "a");
		if (times == NULL) {
			printf("Error opening file");
		exit(1); 
    }
    */
		double *sort; // = (double *)malloc(tries*sizeof(double));
    if(rank==0){
      sort = (double *)malloc(tries*sizeof(double));
		  if(sort==NULL){
			  printf("Memory allocation failed!!\n");
			  free(sort);
			  exit(1);
	    }
    }
 	float** M = NULL;//declareMatrix(n);
	float** T = NULL;//declareMatrix(n);
// 
for(int i=0; i<tries; i++){
	//cache_flush();
	//float** M = NULL;//declareMatrix(n);
	//float** T = NULL;//declareMatrix(n);
   M=declareMatrix(n);
   T=declareMatrix(n);
   //initializeMatrix(M, n); 
	if (rank == 0) {
       // M = declareMatrix(n); 
       initializeMatrix(M, n); 
		   // T = declareMatrix(n);
        s = MPI_Wtime();
    }
  //if(mpicheckSym(M, n ,rank, size)==0){
    mpiMatrixTranspose(M, T, n, rank, size);
 // }
 // else{printf("\nERRORRE\n\n");}
   
	if(rank==0){
		e = MPI_Wtime();
    //if(n<=32){
      //checkCorrectness(M,T,n);
    //}
    
		const char *times_table="times_table.csv";
		FILE *times = fopen(times_table, "a");
		if (times == NULL) {
			printf("Error opening file");
		exit(1); 
    }
    
		time_elapsed = (e-s);
		sort[i]=time_elapsed;
		fprintf(times, "MPI,%d,%d,%.9f\n",n, size, time_elapsed);
    //printf("\nTime elapsed: %.9f\n", time_elapsed);		
	  //freeMatrix(M);
	  //freeMatrix(T);
     fclose(times);
    }
    freeMatrix(M);
    freeMatrix(T);
}
// --------------------------
//sort and average
	//fclose(times);
/*	if(rank==0){
 	  bubbleSort(sort, tries);
    time_avg = timeAverage(sort, tries);
  }

  const char *times_average="times_average.csv";
	FILE *avg = fopen(times_average, "a+");
	if (avg == NULL) {
		printf("Error opening file");
		return 1;
  }*/

if(rank==0){
    bubbleSort(sort, tries);
    //for(int i=0; i<tries; i++){printf("%.9f ", sort[i]);}
    time_avg = timeAverage(sort, tries);
    //printf("time average: %.9f\n", time_avg);
  	
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
      //printf("..");
				gettime=atof(token);
			}
			token = strtok(NULL, ",");
		}
	}
	double speed=0.0;
	double eff=0.0;
 //printf("time: %.9f ", gettime);
	//speed = speedup(gettime, time_avg);  //computes the speedup
	//eff = efficiency(speed, size);  //computes the efficiency
	//printf("\nGettime:%.9f\nAvgTime:%.9f\nSpeedup:%.9f\nEfficiency:%.9f", gettime, time_avg, speed, eff);

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
   //printf("end");
//	  free(sort); 
    //fclose(avg);
    MPI_Finalize();
    return 0;
}
