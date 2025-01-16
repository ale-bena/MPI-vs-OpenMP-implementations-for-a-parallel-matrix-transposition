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
        printf("allocation error");
        exit(1);
    }
    M[0] = (float*)malloc(n * n * sizeof(float));
    if(M[0] == NULL){
        fprintf(stderr, "Error in matrix allocation(2)\n");
        printf("allocation error");
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

int mpicheckSym(float** M, int n, int rank, int size){
    /*float* local_matrix = (float*)malloc(n * n * sizeof(float));
    if (rank == 0) {
        memcpy(local_matrix, M[0], n * n * sizeof(float)); 
    }*/
    //MPI_Bcast(local_matrix, n * n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    if(M==NULL || M[0]==NULL){printf("null in checksym");}
    MPI_Bcast(M?M[0]:NULL, n * n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    int rows_per_process = n/size; 
    int start_row = rank*rows_per_process; 
    int end_row = start_row+rows_per_process;
    if(rank==size-1){end_row=n;}
    // Local check for symmetry
    int local_symmetry = 1;
    for (int i = start_row; (i < end_row) && local_symmetry; i++) {
        for (int j=0; (j<i) && local_symmetry; j++) {
            //if (M[i * n + j] != M[j * n + i]) {
            if (M[i][j] != M[j][i]) {
                local_symmetry = 0;
                //break;
            }
        }
    }
    
    // Combine results from all processes
    int global_symmetry;
    MPI_Reduce(&local_symmetry, &global_symmetry, 1, MPI_INT, MPI_LAND, 0, MPI_COMM_WORLD);
//provare reduce all
    // Broadcast the result back to all processes
    MPI_Bcast(&global_symmetry, 1, MPI_INT, 0, MPI_COMM_WORLD);
    //free(local_matrix);
    return global_symmetry;
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

void localTranspose(float* local_matrix, float* local_transposed, int rows_p, int n) {
    for (int i = 0; i < rows_p; i++) {
        for (int j = 0; j < n; j++) {
            local_transposed[j * rows_p + i] = local_matrix[i * n + j];
        }
    }
}

// Function to handle MPI setup, scatter, gather, and transposition
void mpiMatrixTranspose(float** M, float** T, int n, int rank, int size){
   // float** matrix = NULL;  // Original matrix
   // float** result = NULL;  // Transposed matrix

     // Number of columns per process

   /* print the original matrix
   if (rank == 0) {
        printf("Original matrix:\n");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                printf("%.2f ", M[i][j]);
            }
            printf("\n");
        }
    }*/
    if (mpicheckSym(M, n, rank, size) == 1) {
        if (rank == 0) {
            printf("Matrix is symmetric. No transpose needed.\n");
        }
        return;
    }
    
    
 // if(mpicheckSym(M, n, rank, size)==0){
    int rows_p = n / size;
    if (rank == size - 1) {
      rows_p += n % size;  // Handle remainder
    }
    // Allocate local buffers for scattering and transposing
    float* local_matrix = (float*)malloc(rows_p * n * sizeof(float));
    float* local_transposed = (float*)malloc(rows_p * n * sizeof(float));
    if (local_matrix == NULL || local_transposed == NULL) {
    printf("Failed to allocate local buffers on rank %d.\n", rank);
    //MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Scatter columns from the original matrix to all processes
    MPI_Scatter(M?M[0]:NULL, rows_p * n, MPI_FLOAT, local_matrix, rows_p * n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    //printf("scatter recived %d ", rank );
/*    // Print the scattered rows for debugging
    for (int rank_id = 0; rank_id < size; rank_id++) {
        if (rank == rank_id) {
            printf("\n------ RANK %d RECEIVED ROWS -----\n", rank);
            for (int i = 0; i < rows_p; i++) {
                for (int j = 0; j < n; j++) {
                    printf("%.2f ", local_matrix[i * n + j]);
                }
                printf("\n");
            }
        }
        MPI_Barrier(MPI_COMM_WORLD); // Synchronize for ordered printing
        
    }
 */ 
 
    for (int i = 0; i < rows_p; i++) {
        for (int j = 0; j < n; j++) {
            local_transposed[j * rows_p + i] = local_matrix[i * n + j];
        }
    }
    //printf("transpose done %d ", rank);
    //localTranspose(local_matrix, local_transposed, rows_p, n);
    /*
 
    for (int i = 0; i < rows_p; i++) {
        for (int j = 0; j < n; j++) {
            local_transposed[j * rows_p + i] = local_matrix[i * n + j];
        }
    }
    */

/*  MPI_Barrier(MPI_COMM_WORLD);
    
    // Print the scattered columns for debugging
    for (int rank_id = 0; rank_id < size; rank_id++) {
        if (rank == rank_id) {
            printf("\n------ RANK %d RECEIVED COLS -----\n", rank);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < rows_p; j++) {
                    printf("%.2f ", local_transposed[i * rows_p + j]);
                }
                printf("\n");
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (rank == 0) {
        result = declareMatrix(n);  //declared in the main
    }
*/
    //MPI_Barrier(MPI_COMM_WORLD);
// Gather transposed columns back to the master process
    MPI_Gather(local_transposed, rows_p * n, MPI_FLOAT,
               T ? T[0] : NULL, rows_p * n, MPI_FLOAT,
               0, MPI_COMM_WORLD);
// Reorganize the gathered data on the root process
   if (rank == 0) {
     //printf(" reordering ");
        float* temp = (float*)malloc(n * n * sizeof(float));
        for (int r = 0; r < size; r++) {
            //int start = r * rows_p * n;
            for(int i=0; i<n; i++){
               for(int j=0; j<rows_p; j++){
                   temp[i * n + j + r*rows_p] = T[0][ r * rows_p * n + i*rows_p + j];
               }              
            }
        }
        //printf("copying");
// Copy the reorganized data back into the result matrix
        memcpy(T[0], temp, n * n * sizeof(float));
        //printf("done");
        /* Print the transposed matrix
        printf("\nTransposed matrix:\n");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                printf("%.2f ", T[i][j]);
            }
            printf("\n");
        }*/

        // Free the matrices
       // freeMatrix(matrix);
       // freeMatrix(result);
       free(temp);
    }
  //printf("freeing");
    // Free local buffers

    free(local_matrix);
    free(local_transposed);
    //printf("finished");
    
  //} //if checkSym parenthesis
  //else{printf("ERRORRRRRRR!");}
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
        array[i] = 1; 
    }

    free(array);
}

