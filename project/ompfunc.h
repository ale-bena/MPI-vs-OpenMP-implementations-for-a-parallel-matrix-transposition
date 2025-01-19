#ifndef OMPFUNC_H_
#define OMPFUNC_H_

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>


int isPowerOfTwo(int n); //verifies if the number is a power of 2

float random_float(int min, int max); //generates random float numbers

float** declareMatrix(int n); //declares the matrix

void initializeMatrix(float** M, int n); //initializes the matrix with random numbers

int checkSymOMP(float** M, int n); //checks if the matrix is symmetric

void matTransposeOMP(float** M, float** T, int n); //transposes the matrix

void printMatrix(float **M, int n); //prints the matrix, to see if the program works

//void freeMemory(float** M, int n); //frees the allocated memory

void freeMemory(float** M);

void bubbleSort(double sort[], int n);  //bubble sort function for a double array

double timeAverage(double sort[], int n); //computes the average of the 40% of the elements of the array

void checkCorrectness(float **M, float** T, int n);  //checks if the transposition has been executed correctly

double speedup(double timeSeq, double timeAvg);  //computes the speedup

double efficiency(double speedup, int nthreads);  //computes the efficiency

void cache_flush(void);  //frees the cache


#endif
