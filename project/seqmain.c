#include "seqfunc.h"

int main(int argc, char** argv){
	if(argc<2){
		printf("Insert a number!!");
		exit(1);
	}
	int n=atoi(argv[1]);
	if(isPowerOfTwo(n)==0){
		exit(1);
	}
	
	int tries=atoi(argv[2]);
	if(tries<1){
		exit(1);
	}
	
	srand(25); //random fixed seed
	
	double s, e, time_elapsed; //omp time measurement

	const char *times_table="times_table.csv";  //file for time values
	
	FILE *times = fopen(times_table, "a");
	if (times == NULL) {
		perror("Error opening file");
		return 1;
	}
	
	fseek(times, 0, SEEK_END);
	if(ftell(times)==0){  //inserting the header
		fprintf(times, "%-10s,%-10s,%-10s,%-10s\n", "Type", "Dim", "Threads", "Time");
	}

	double *sort = (double *)malloc(tries*sizeof(double)); //sort array initialization
	if(sort==NULL){
		printf("Memory allocation failed!!\n");
		free(sort);
		exit(1);
	}
  float** M=NULL;
  float** T=NULL;
	for(int i=0; i<tries; i++) {
		cache_flush();
		M=declareMatrix(n);
		initializeMatrix(M, n);
		T=declareMatrix(n);
		s=omp_get_wtime();
		matTranspose(M,T,n);
		e=omp_get_wtime();
    if(n<=32){
      checkCorrectness(M,T,n);
    }
		time_elapsed = (e-s);
		sort[i]=time_elapsed;
		fprintf(times, "Sequential,%d,1,%.9f\n", n, time_elapsed);		
		freeMemory(M);
		freeMemory(T);
    M=NULL; T=NULL;
	}
	fclose(times);
	
	bubbleSort(sort, tries);
	
	double time_avg = timeAverage(sort, tries); 
	
	const char *times_average="times_average.csv";
	
	//const char *header2 = "Type, Dim, Threads, Time, Speedup, Efficiency";
	FILE *avg = fopen(times_average, "a");
	if (avg == NULL) {
		perror("Error opening file");
		return 1;
	}
	
	fseek(avg, 0, SEEK_END);
	if(ftell(avg)==0){
		//fprintf(avg, "%s\n", header2);
fprintf(avg, "%-10s,%-10s,%-10s,%-10s,%-10s,%-10s\n", "Type", "Dim", "Threads", "Average time", "Speedup", "Efficiency");
//fprintf(avg, "%-10s%-10s%-10s%-10s,%-10s,%-10s\n", "Type,", " Dim,", " Threads,", "Average time", "Speedup", "Efficiency");
	}
	fprintf(avg, "Sequential,%d,1,%.9f,1,1\n", n, time_avg);
	
	fclose(avg);
	free(sort); 
    
    return 0;
}
