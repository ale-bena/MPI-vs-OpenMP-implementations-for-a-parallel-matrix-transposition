#ifndef MPIFUNC_H_
#define MPIFUNC_H_

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>

int isPowerOfTwo(int n); //verifies if the number is a power of 2

float random_float(int min, int max); //generates random float numbers

float** declareMatrix(int n); //declares the matrix

void initializeMatrix(float** M, int n); //initializes the matrix with random numbers

void printMatrix(float **M, int n); //prints the matrix, to see if the program works

void bubbleSort(double sort[], int n);  //bubble sort function for a double array

double timeAverage(double sort[], int n); //computes the average of the 40% of the elements of the array

void checkCorrectness(float **M, float** T, int n);  //checks if the transposition has been executed correctly

double speedup(double timeSeq, double timeAvg);  //computes the speedup

double efficiency(double speedup, int nthreads);  //computed the efficiency

void cache_flush();  //frees the cache

void matTransposeMPI(float** M, float**T, int n, int rank, int size);  //transpose function

//void localTranspose(float* local_columns, float* local_transposed, int cols_per_process, int n);

int checkSymMPI(float** M, int n, int rank, int size);  //checks if the matrix is symmetric

void freeMatrix(float** M);  //frees allocated memory for the matrix



#endif
