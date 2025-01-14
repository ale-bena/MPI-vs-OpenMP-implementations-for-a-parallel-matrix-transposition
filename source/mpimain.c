#include "mpifunc.h"
#define MAX 128

int main(int argc, char** argv){
    int rank, size;
    int n = 16; // Matrix size (n x n)
    srand(25);
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Ensure the number of processes divides n evenly
    if (n % size != 0 && rank == 0) {
        fprintf(stderr, "Matrix size %d must be divisible by the number of processes %d.\n", n, size);
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    // Perform matrix transposition with MPI
    mpiMatrixTranspose(n, rank, size);

    MPI_Finalize();

/*	if(argc<2){
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
	double s,e, time_elapsed;
	const char *times_table="times_table.csv";
	const char *header1 = "Type, Dim, Threads, Time";
	FILE *times = fopen(times_table, "a");
	if (times == NULL) {
		perror("Error opening file");
    exit(1); 
	}

//SORT ARRAY INITIALIZATION
	double *sort = (double *)malloc(tries*sizeof(double));
	if(sort==NULL){
		printf("Memory allocation failed!!\n");
		free(sort);
		exit(1);
	}
 
   switch(n){
     case 16: blocksize=8; break;
     case 32: blocksize=16; break;
     case 64: blocksize=32; break;
     default: blocksize=64;
   }
 
//EXECUTION
	for (int i=0; i<tries; i++){
	  cache_flush();
		float** M=declareMatrix(n);
		initializeMatrix(M, n);
		float** T=declareMatrix(n);
		s=omp_get_wtime();
		matTransposeOMP(M,T,n,blocksize);
		e=omp_get_wtime();
		time_elapsed = (e-s);
		sort[i]=time_elapsed;
		fprintf(times, "BlockOmp,%d,%d,%.9f\n",n, n_threads, time_elapsed);		
		freeMemory(M, n);
		freeMemory(T,n);
	} 
  fclose(times);
 	bubbleSort(sort, tries);
 	double time_avg = timeAverage(sort, tries); 
  const char *times_average="times_average.csv";
	const char *header2 = "Type, Dim, Threads, Time, Speedup, Efficiency";
	FILE *avg = fopen(times_average, "a+");
	if (avg == NULL) {
		perror("Error opening file");
		return 1;
	}

	char buffer[MAX];
	char string[] = "Sequential";
	double gettime = 0.0;
	fgets(buffer, MAX, times);
	while(fgets(buffer, MAX, times)){
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
	speed = speedup(gettime, time_avg);  //computes the speedup
	eff = efficiency(speed, n_threads);  //computes the efficiency
	//printf("\nGettime:%.9f\nAvgTime:%.9f\nSpeedup:%.9f\nEfficiency:%.9f", gettime, time_avg, speed, eff);
	if(gettime==0.0){ //handles the case where the sequential program is not run before this program
     fprintf(avg, "BlockOmp%d,%d,%d,%.9f,-,-\n", blocksize, n, n_threads, time_avg);
   }
   else{
	  fprintf(avg, "BlockOmp,%d,%d,%.9f,%.9f,%.9f\n", n, n_threads, time_avg, speed, eff);
	}
	fclose(avg);
	free(sort); 
 */
 
    return 0;
}
