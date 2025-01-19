#include "mpifunc.h"
#define min(a, b) ((a) < (b) ? (a) : (b))

int isPowerOfTwo(int n) {
    return (n > 0) && ((n & (n - 1)) == 0);
}

float random_float(int min, int max){
	float n;
	int tmp = 0;
	tmp = rand() % ((max-min+1)+min);
	n = (float)tmp/100;
	return n;
}

float** declareMatrix(int n){
    float** M = (float**)malloc(n * sizeof(float*));
    if(M == NULL){
        printf("Error in matrix allocation\n");
        printf("allocation error");
        exit(1);
    }
    M[0] = (float*)malloc(n * n * sizeof(float));
    if(M[0] == NULL){
        printf("Error in matrix allocation(2)\n");
        free(M);
        exit(1);
    }
    for(int i=1; i<n; i++){
        M[i] = M[0] + i * n;
    }
    return M;
}


void initializeMatrix(float** M, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            M[i][j] = random_float(0,9999);
        }
    }
}


int checkSymMPI(float** M, int n, int rank, int size){
    MPI_Bcast(M[0], n * n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    int rows_per_process = n / size;
    int start_row = rank * rows_per_process;
    int end_row = (rank == size - 1) ? n : start_row + rows_per_process;
    // Local check for symmetry
    int local_symmetry = 1;
    for (int i = start_row; (i < end_row) && local_symmetry; i++) {
        for (int j=0; (j<i) && local_symmetry; j++) {
            if (M[i][j] != M[j][i]) {
                local_symmetry = 0;
            }
        }
    }
    int global_symmetry;
    MPI_Allreduce(&local_symmetry, &global_symmetry, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    return global_symmetry;
}



/*
void localTranspose(float* local_matrix, float* local_transposed, int rows_p, int n) {
    for (int i = 0; i < rows_p; i++) {
        for (int j = 0; j < n; j++) {
            local_transposed[j * rows_p + i] = local_matrix[i * n + j];
        }
    }
}*/



void matTransposeMPI(float** M, float** T, int n, int rank, int size){
  if(mpicheckSym(M, n, rank, size)==0){
    int rows_p = n / size;
    if (rank == size - 1) {
      rows_p += n % size;  // Handle remainder
    }
    float* temp=NULL;
    if(rank==0){temp = (float*)malloc(n * n * sizeof(float));}
    /*float* local_matrix = (float*)malloc(rows_p * n * sizeof(float));
    float* local_transposed = (float*)malloc(rows_p * n * sizeof(float));
    if (local_matrix == NULL || local_transposed == NULL) {
      printf("Failed to allocate local buffers on rank %d.\n", rank);
    }*/
    MPI_Scatter(M?M[0]:NULL, rows_p * n, MPI_FLOAT, T[0], rows_p * n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    for (int i = 0; i < rows_p; i++) {
        for (int j = 0; j < n; j++) {
            //local_transposed[j * rows_p + i] = local_matrix[i * n + j];
            M[0][j * rows_p + i] = T[0][i * n + j];
        }
    }
    MPI_Gather(M[0], rows_p * n, MPI_FLOAT,temp, rows_p * n, MPI_FLOAT,0, MPI_COMM_WORLD);      
// Reorganize the gathered data on the root process
   if (rank == 0) {
        for (int r = 0; r < size; r++) {
            for(int i=0; i<n; i++){
               for(int j=0; j<rows_p; j++){
                   T[0][i * n + j + r*rows_p] = temp[ r * rows_p * n + i*rows_p + j];
               }              
            }
        }
       free(temp);
    }
   // free(local_matrix);
   // free(local_transposed);    
  }
  else{
    printf("MATRIX IS SYMMETRIC");
  }
}

void printMatrix(float **M, int n){
 for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
        printf("%.2f  ",M[i][j]);
    }
    printf("\n");
 }
}


void freeMatrix(float** M) {
    free(M[0]);
    free(M); 
}

void bubbleSort(double sort[], int n) { 
    for (int i = 0; i < n - 1; i++) { 
        for (int j = 0; j < n - i - 1; j++) { 
            if (sort[j] > sort[j + 1]) { 
                double temp = sort[j]; 
                sort[j] = sort[j + 1]; 
                sort[j + 1] = temp; 
            } 
        } 
    } 
} 
 
 
double timeAverage(double sort[], int n) { 
    int range = (int)(n * 0.4); 
    if (range == 0) 
    	return 0.0; 
 
    int starti = (int)((n - range) / 2);  
    int endi = starti + range;  
 
    double sum = 0.0; 
    for (int i = starti; i < endi; i++) { 
    	//printf("sort%d: %.9f", i, sort[i]);
        sum += sort[i]; 
    } 
    return sum / range; 
} 

void checkCorrectness(float **M, float** T, int n){ 
    int tmp = n/4;
    int i=1;
    for(int j=1; j<n; j+=tmp){
      if(M[i][j]!=T[j][i]){
        printf("\n\nError in the transposition!!\n\n");
      }
    }
}

double speedup(double timeSeq, double timeAvg){
	double tmp = 0.0;
	tmp = timeSeq/timeAvg;
	//printf("\nTMP1:%.9f", tmp);
	return tmp; 
}

double efficiency(double speedup, int nthreads){
	double tmp = 0.0;
	tmp = (speedup/nthreads)*100;
	//printf("\nTMP2:%.9f\n", tmp);
	return tmp;
}

void cache_flush() {
 	size_t cache_size = 1088*1024;
	char *array = (char*)malloc(cache_size); 
    if (array == NULL) {
        printf("Memory allocation failed\n");
        exit(1);
    }

    for (size_t i = 0; i < cache_size; i++) {
        array[i] = 1; 
    }

    free(array);
}

