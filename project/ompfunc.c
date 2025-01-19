#include "ompfunc.h"

int isPowerOfTwo(int n) {
    return (n > 0) && ((n & (n - 1)) == 0); 
}

float random_float(int min, int max){
	float n;
	int tmp=0;
	tmp = rand() % ((max-min+1)+min);
	n=(float)tmp/100;
	return n;
}

/*
float** declareMatrix(int n){
	float** M = (float**)malloc(n * sizeof(float*));
	if(M==NULL){
		exit(1);
	}
	for (int i = 0; i < n; ++i) {
        	M[i] = (float*)malloc(n * sizeof(float));
        	if(M[i]==NULL){
        		free(M);
        		exit(1);
        	}
    	}
    	return M;
}
*/
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
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            M[i][j] = random_float(0,9999);
        }
    }
}


/***************************************
*
* checkSymOMP() - checks if a given n x n 
* matrix is symmetric using OpenMP for parallelism
* input: float** M - pointer to the matrix
*        int n - the size of the matrix 
*        (number of rows and columns)
* output: 1 (true) if the matrix is symmetric,
*         0 (false) otherwise
* notes: uses OpenMP parallelization with a logical 
*        AND reduction to efficiently combine results 
*        from multiple threads. Symmetry is determined 
*        by comparing M[i][j] with M[j][i].
*
***************************************/
int checkSymOMP(float** M, int n) {
    int tmp = 1;
#pragma omp parallel for reduction(&&:tmp)
    for (int i = 1; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            if (M[i][j] != M[j][i]) {
                tmp = 0; 
            }
        }
    }
    return tmp;
}


/*void matTransposeOMP(float** M, float** T, int n){
    if(checkSym(M,n)==0){
#pragma omp for collapse(2)
        	for (int i = 0; i < n; ++i) {
            		for (int j = 0; j < n; ++j) {
                		T[i][j] = M[j][i];
        	}
        }
    }
}*/

/***************************************
*
* matTransposeOMP() - transposes an n x n 
* matrix using OpenMP for parallelism with 
* prefetching for optimization
* input: float** M - pointer to the original matrix
*        float** T - pointer to the transposed matrix
*        int n - the size of the matrix 
*        (number of rows and columns)
* output: none (matrix T is modified in place)
* notes: checks if the matrix M is symmetric using 
*        checkSymOMP(). If not symmetric, performs 
*        parallel transposition using OpenMP. Prefetching 
*        is applied to improve memory access efficiency.
*
***************************************/
void matTransposeOMP(float** M, float** T, int n) {
    if (checkSymOMP(M, n) == 0) {
#pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (j + 1 < n){
                    __builtin_prefetch(&M[j + 1][i], 0, 1); // Prefetch for read
                    __builtin_prefetch(&T[i][j + 1], 1, 1); // Prefetch for write
                }
                T[i][j] = M[j][i];
            }
        }
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
/*
void freeMemory(float** M, int n){
	for(int i=0; i<n; i++){
		free(M[i]);
	}
	free(M);
}*/

void freeMemory(float** M) {
    free(M[0]);
    free(M); 
}


// Bubble Sort Function 
void bubbleSort(double sort[], int n) { 
    for (int i = 0; i < n - 1; i++) { 
        for (int j = 0; j < n - i - 1; j++) { 
            if (sort[j] > sort[j + 1]) { 
                // Swap sort[j] and sort[j+1] 
                double temp = sort[j]; 
                sort[j] = sort[j + 1]; 
                sort[j + 1] = temp; 
            } 
        } 
    } 
} 
 
// Central time average function 
double timeAverage(double sort[], int n) { 
    int range = (int)(n * 0.4); // Calculate 40% of the total count!!!cambiare il .4 per cambiare tutta la funzione
    if (range == 0) 
    	return 0.0; // Handle edge case where rangeCount is 0 
 
    int starti = (n - range) / 2; // Start of the central range 
    int endi = starti + range; // End of the central range 
 
    double sum = 0.0; 
    for (int i = starti; i < endi; i++) { 
    	//printf("sort%d: %.9f", i, sort[i]);
        sum += sort[i]; 
    } 
    return sum / range; // Calculate the average  guardare che non ci siano problemi di ritorno perche faccio double/int
} 

// Check correctness function
void checkCorrectness(float **M, float** T, int n){ //mettere controllo per n<=32 fuori dalla funzione 
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

// Cache reset function
/*void cache_flush(){
	size_t l1_cache = 32*1024;
	size_t l2_cache = 1024*1024;
	char *array_l1_1 = (char*)malloc(l1_cache); 
	char *array_l1_2 = (char*)malloc(l1_cache); 
	char *array_l2 = (char*)malloc(l2_cache); 
	if(array_l1_1==NULL || array_l1_2==NULL || array_l2==NULL){
		printf("Memory allocation failed\n");
		free(array_l1_1);
		free(array_l1_2);
		free(array_l2);
		exit(1);
	}
  for(size_t i=0; )
	memset(array_l1_1, 0, l1_cache);
	memset(array_l1_2, 0, l1_cache);
	memset(array_l2, 0, l2_cache);
	free(array_l1_1);
	free(array_l1_2);
	free(array_l2);
}*/

void cache_flush() {
 	  size_t cache_size = 1088*1024;
	  char *array = (char*)malloc(cache_size); 
    if (array == NULL) {
        printf("Memory allocation failed\n");
        exit(1);
    }

    // Populate the array with 0, touching each cache line
    for (size_t i = 0; i < cache_size; i++) {
        array[i] = 0; 
    }

    free(array);
}



