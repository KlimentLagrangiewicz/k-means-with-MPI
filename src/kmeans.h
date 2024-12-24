#ifndef KMEANS_H_
#define KMEANS_H_

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "mpi.h"

double getEvDist(const double *x1, const double *x2, int m);
void autoscaling(double* const x, const int n, const int m);
void blockFunction1(const double* const x, double* const Ex, double* const Exx, const int m, const int perProc);
void blockFunction2(double* const x, const double* const Ex, const double* const Exx, const int m, const int perProc);
void MPI_Autoscaling(double* const x, const int n, const int m);
int getCluster(const double* const x, const double* const c, const int m, int k);
char constr(const int *y, const int val, int s);
void detCores(const double* const x, double* const c, const int n, const int m, const int k);
void detStartPartition(const double* const x, const double* const c, int* const y, int* const nums, int n, const int m, const int k);
void MPI_Detstartpartition(const double* const x, const double* const c, int* const y, int* const nums, const int n, const int m, const int k);
void sumCores(const double* const x, double* const c, const int* const y, const int n, const int m);
void MPI_Sumcores(const double* const x, double* const c, const int* const y, const int n, const int m, const int k);
void coresDiv(double* const c, const int* const nums, const int m, const int k);
void MPI_Coresdiv(double* const c, const int* const nums, const int m, const int k);
char checkPartition(const double* const x, const double* const c, int* const y, int* const nums, const int n, const int m, const int k);
char MPI_Checkpartition(const double* const x, const double* const c, int* const y, int* const nums, const int n, const int m, const int k);
char MPI_Cyclicrecalc(const double* const x, double* const c, int* const y, int* const nums, const int n, const int m, const int k);
void MPI_Kmeans(const double* const X, int* const y, const int n, const int m, const int k);

#endif
