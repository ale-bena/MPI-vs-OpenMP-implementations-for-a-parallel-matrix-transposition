#include "mpifunc.h"
#define min(a, b) ((a) < (b) ? (a) : (b))

int isPowerOfTwo(int n) {
    return (n > 0) && ((n & (n - 1)) == 0); //uses the binary 
}

float random_float(int min, int max){
	float n;
	int tmp = 0;
	tmp = rand() % ((max-min+1)+min);
	n = (float)tmp/100;
	return n;
}
/*
float** declareMatrix(int n){ //meglio allocarla in modo continuo
	float** M = (float**)malloc(n * sizeof(float*));
	if(M==NULL){
		exit(1);
	}
	for (int i = 0; i < n; i++) {
        	M[i] = (float*)malloc(n * sizeof(float));
        	if(M[i]==NULL){ //non sto deallocando le righe già allocate
        		free(M);
        		exit(1);
        	}
    	}
    	return M;
}
*/
float** declareMatrix(int n){
    if(n <= 0) {
        return NULL;
    }
    float** M = (float**)malloc(n * sizeof(float*));
    if(M == NULL){
        fprintf(stderr, "Error in matrix allocation\n");
        exit(1);
    }
    M[0] = (float*)malloc(n * n * sizeof(float));
    if(M[0] == NULL){
        fprintf(stderr, "Error in matrix allocation(2)\n");
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

/***************************************
*
* checkSymOMP() - checks if a given n x n 
* matrix is symmetric using OpenMP for parallelism
* input: float** M - pointer to the matrix
*        int n - the size of the matrix 
*        (number of rows and columns)
* output: 1 (true) if the matrix is symmetric,
*         0 (false) otherwise
* notes: utilizes OpenMP parallelization with a 
*        reduction operation to combine results from
*        multiple threads. Symmetry is determined by 
*        comparing M[i][j] with M[j][i].
*
***************************************/
int checkSymOMP(float** M, int n) {
    int tmp = 1;
#pragma omp parallel for reduction(&&:tmp)
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < i; j++) {
            if (M[i][j] != M[j][i]) {
                tmp = 0; 
            }
        }
    }
    return tmp;
}

// Function to check if a matrix is symmetric
int checkSym(float** M, int n) {
    for (int i = 1; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            if (M[i][j] != M[j][i]) {
                return 0;
            }
        }
    }
    return 1;
}



/***************************************
*
* matTransposeOMP() - transposes an n x n 
* matrix using OpenMP for parallelism and 
* block-based optimization
* input: float** M - pointer to the original matrix
*        float** T - pointer to the transposed matrix
*        int n - the size of the matrix 
*        (number of rows and columns)
*        int blockSize - size of the blocks used for 
*        block-based transposition
* output: none (matrix T is modified in place)
* notes: checks if the matrix M is symmetric using 
*        checkSymOMP(). If not symmetric, performs 
*        block-based transposition with parallelization 
*        to optimize performance on large matrices.
*
***************************************/
void matTransposeOMP(float** M, float** T, int n, int blockSize) {
    if(checkSymOMP(M,n)==0){
#pragma omp parallel for collapse(2) shared(M, T)
    for (int i = 0; i < n; i += blockSize) {
        for (int j = 0; j < n; j += blockSize) {
            int check=min(n, (i+blockSize));
            int check2=min(n, (j+blockSize));
            for (int ii = i; ii < check; ii++) {
                for (int jj = j; jj < check2; jj++) {
                    T[ii][jj] = M[jj][ii];
                }
            }
        }
    }
    }
}

// Function to transpose a local chunk of columns with symmetry check
void matTranspose(float** M, float** T, int n){
    if (checkSym(M, n)==0) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                T[i][j] = M[j][i];
            }
        }
    }
}

void localTranspose(float* local_columns, float* local_transposed, int cols_per_process, int n) {
    for (int i = 0; i < cols_per_process; i++) {
        for (int j = 0; j < n; j++) {
            local_transposed[j * cols_per_process + i] = local_columns[i * n + j];
        }
    }
}

// Function to handle MPI setup, scatter, gather, and transposition
void mpiMatrixTranspose(int n, int rank, int size){
    float** matrix = NULL;  // Original matrix
    float** result = NULL;  // Transposed matrix

    int cols_per_process = n / size; // Number of columns per process

    // Master process initializes the matrix
    if (rank == 0) {
        matrix = declareMatrix(n);
        initializeMatrix(matrix, n);

        /*for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                matrix[j][i] = i * n + j; // Fill with some values
            }
        }*/

        printf("Original matrix:\n");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                printf("%.2f ", matrix[i][j]); // Access by rows for display
            }
            printf("\n");
        }
    }

    // Allocate local buffers for scattering and transposing
    float* local_columns = (float*)malloc(cols_per_process * n * sizeof(float));
    float* local_transposed = (float*)malloc(cols_per_process * n * sizeof(float));

    // Scatter columns from the original matrix to all processes
    MPI_Scatter(matrix ? matrix[0] : NULL, cols_per_process * n, MPI_FLOAT,
                local_columns, cols_per_process * n, MPI_FLOAT,
                0, MPI_COMM_WORLD);

  /*    // Print the scattered columns for debugging
    for (int rank_id = 0; rank_id < size; rank_id++) {
        if (rank == rank_id) {
            printf("\n------ RANK %d RECEIVED COLS -----\n", rank);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < cols_per_process; j++) {
                    printf("%.2f ", local_columns[i * n + j]);
                }
                printf("\n");
            }
        }
        MPI_Barrier(MPI_COMM_WORLD); // Synchronize for ordered printing
    }*/
    // Print the scattered rows for debugging
    for (int rank_id = 0; rank_id < size; rank_id++) {
        if (rank == rank_id) {
            printf("\n------ RANK %d RECEIVED ROWS -----\n", rank);
            for (int i = 0; i < cols_per_process; i++) {
                for (int j = 0; j < n; j++) {
                    printf("%.2f ", local_columns[i * n + j]);
                }
                printf("\n");
            }
        }
        MPI_Barrier(MPI_COMM_WORLD); // Synchronize for ordered printing
        
    }
  
    localTranspose(local_columns, local_transposed, cols_per_process, n);

    MPI_Barrier(MPI_COMM_WORLD);
    
    // Print the scattered columns for debugging
    for (int rank_id = 0; rank_id < size; rank_id++) {
        if (rank == rank_id) {
            printf("\n------ RANK %d RECEIVED COLS -----\n", rank);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < cols_per_process; j++) {
                    printf("%.2f ", local_transposed[i * cols_per_process + j]);
                }
                printf("\n");
            }
        }
        MPI_Barrier(MPI_COMM_WORLD); // Synchronize for ordered printing
    }

    if (rank == 0) {
        result = declareMatrix(n); 
    }

    // Gather transposed columns back to the master process
    MPI_Gather(local_transposed, cols_per_process * n, MPI_FLOAT,
               result ? result[0] : NULL, cols_per_process * n, MPI_FLOAT,
               0, MPI_COMM_WORLD);


 // Reorganize the gathered data on the root process
   if (rank == 0) {
        float* temp = (float*)malloc(n * n * sizeof(float));
        for (int r = 0; r < size; r++) {
            //int start_col = r * cols_per_process;
            int start = r * cols_per_process * n;
            for(int i=0; i<n; i++){
               for(int j=0; j<cols_per_process; j++){
                   temp[i * n + j + r*cols_per_process] = result[0][start + i*cols_per_process + j];
               }              
            }
            
            /*for (int i = 0; i < cols_per_process; i++) {
                for (int j = 0; j < n; j++) {
                    temp[j * n + start_col + i] = result[0][r * cols_per_process * n + i * n + j];
                }
            }*/
        }

        // Copy the reorganized data back into the result matrix
        memcpy(result[0], temp, n * n * sizeof(float));
        free(temp);

        // Print the transposed matrix
        printf("\nTransposed matrix:\n");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                printf("%.2f ", result[i][j]);
            }
            printf("\n");
        }

        // Free the matrices
        freeMatrix(matrix);
        freeMatrix(result);
    }

    // Free local buffers
    free(local_columns);
    free(local_transposed);
}


void printMatrix(float **M, int n){
 for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
        printf("%.2f  ",M[i][j]);
    }
    printf("\n");
 }
}


void freeMemory(float** M, int n){ //sistemare per continous allocation
	for(int i=0; i<n; i++){
		free(M[i]);
	}
	free(M);
}

void freeMatrix(float** M) {
    free(M[0]); // Free contiguous block
    free(M);    // Free column pointers
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
 
    int starti = (n - range) / 2;  
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
        array[i] = 0; 
    }

    free(array);
}

