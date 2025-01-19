//SEQUENTIAL CODE

#include "seqfunc.h"

/***************************************
*
* isPowerOfTwo() - checks if a number 
* is a power of two
* input: int n - the number to check
* output: 1 (true) if n is a power of two,
*         0 (false) otherwise
*
***************************************/
int isPowerOfTwo(int n) {
    return (n > 0) && ((n & (n - 1)) == 0); 
}

/***************************************
*
* random_float() - generates a random 
* floating-point number within a specified 
* range
* input: int min - the minimum value (inclusive)
*        int max - the maximum value (inclusive)
* output: float n - a random float between 
*         min and max (scaled by 1/100)
*
***************************************/
float random_float(int min, int max){
	float n;
	int tmp=0;
	tmp = rand() % ((max-min+1)+min);
	n=(float)tmp/100;
	return n;
}

/***************************************
*
* declareMatrix() - dynamically allocates 
* memory for an n x n matrix of floats
* input: int n - the size of the matrix 
*        (number of rows and columns)
* output: float** M - pointer to the 
*         allocated matrix
* notes: exits the program with code 1 if 
*        memory allocation fails
*
***************************************/

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
}*/

/***************************************
*
* initializeMatrix() - populates an n x n 
* matrix with random floating-point values
* input: float** M - pointer to the matrix 
*        int n - the size of the matrix 
*        (number of rows and columns)
* output: none (matrix M is modified in place)
* notes: uses random_float() to generate 
*        values in the range [0, 9999]
*
***************************************/
void initializeMatrix(float** M, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
        	M[i][j] = random_float(0,9999);
        }
    }
}

/***************************************
*
* checkSym() - checks if a given n x n 
* matrix is symmetric
* input: float** M - pointer to the matrix
*        int n - the size of the matrix 
*        (number of rows and columns)
* output: 1 (true) if the matrix is symmetric,
*         0 (false) otherwise
* notes: symmetry is determined by comparing
*        M[i][j] with M[j][i] for all i > j
*
***************************************/
int checkSym(float** M, int n){ 
    for (int i = 1; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            if(M[i][j]!=M[j][i]){
                return 0;
            }
        }
    }
    return 1;
}

/***************************************
*
* matTranspose() - transposes a given n x n 
* matrix if it is not symmetric
* input: float** M - pointer to the original matrix
*        float** T - pointer to the matrix to store 
*        the transposed result
*        int n - the size of the matrix 
*        (number of rows and columns)
* output: none (matrix T is modified in place)
* notes: if the matrix M is symmetric (checked using 
*        checkSym()), the transpose operation is skipped
*
***************************************/
void matTranspose(float** M, float** T, int n){
    if(checkSym(M,n)==0){
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                T[i][j] = M[j][i];
            }
        }
    }
}

// used for debugging
void printMatrix(float **M, int n){
 for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
        printf("%.2f  ",M[i][j]);
    }
    printf("\n");
 }
}

/***************************************
*
* freeMemory() - frees the dynamically 
* allocated memory for an n x n matrix
* input: float** M - pointer to the matrix
*        int n - the size of the matrix 
*        (number of rows and columns)
* output: none (frees memory for matrix M)
* notes: frees each row and then the matrix itself
*
***************************************/
/*void freeMemory(float** M, int n){
	for(int i=0; i<n; i++){
		free(M[i]);
	}
	free(M);
}
*/
void freeMemory(float** M) {
    free(M[0]);
    free(M); 
}

/***************************************
*
* bubbleSort() - sorts an array of doubles 
* in ascending order using the bubble sort algorithm
* input: double sort[] - the array to be sorted
*        int n - the number of elements in the array
* output: none (the array is sorted in place)
* notes: repeatedly swaps adjacent elements if they 
*        are in the wrong order, iterating through the array
*
***************************************/
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
 
/***************************************
*
* timeAverage() - calculates the average of the 
* central 40% values in a sorted array of doubles
* input: double sort[] - the sorted array
*        int n - the number of elements in the array
* output: double - the average of the central 40% of the array
* notes: the function calculates 40% of the array size, 
*        determines the central range, and computes the average 
*        of the values within that range. If the range is 0, 
*        returns 0.0.
*
***************************************/
double timeAverage(double sort[], int n) { 
    int range = (int)(n * 0.4); 
    if (range == 0) 
    	return 0.0;  
 
    int starti = (n - range) / 2; 
    int endi = starti + range;  
 
    double sum = 0.0; 
    for (int i = starti; i < endi; i++) { 
    	//printf("sort%d: %.9f", i, sort[i]);
        sum += sort[i]; 
    } 
    return sum / range;   
} 

/***************************************
*
* checkCorrectness() - checks the correctness of the 
* transposition of a matrix by comparing selected elements
* between the original matrix M and its transposed matrix T
* input: float** M - pointer to the original matrix
*        float** T - pointer to the transposed matrix
*        int n - the size of the matrix (number of rows and columns)
* output: none (prints an error message if the transposition is incorrect)
* notes: compares elements at regular intervals to detect errors in the 
*        transposition. If a mismatch is found, it prints an error message.
*
***************************************/
void checkCorrectness(float **M, float** T, int n){  
    int tmp = n/4;
    int i=1;
    for(int j=1; j<n; j+=tmp){
      if(M[i][j]!=T[j][i]){
        printf("\n\nError in the transposition!!\n\n");
      }
    }
}
 
// Computes the speedup 
double speedup(double timeSeq, double timeAvg){
	double tmp = 0.0;
	tmp = timeSeq/timeAvg;
	return tmp; 
}

// Computes the efficiency
double efficiency(double speedup, int nthreads){
	double tmp = 0.0;
	tmp = (speedup/nthreads)*100;
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

    for (size_t i = 0; i < cache_size; i++) {
        array[i] = 0; 
    }

    free(array);
}


