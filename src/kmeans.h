#ifndef KMEANS_H_
#define KMEANS_H_

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include "mpi.h"

double getEvDist(const double *x1, const double *x2, const int m);
void autoscaling(double* const x, const int n, const int m);
void blockFunction1(const double* const x, double* const Ex, double* const Exx, const int m, const int perProc, const int rowID);
void blockFunction2(double* const x, const double* const Ex, const double* const Exx, const int m, const int perProc, const int rowID);
void MPI_Autoscaling(double* const x, const int n, const int m);
int getCluster(const double *x, const double *c, const int m, const int k);
void detCores(const double *x, double *c, const int *sn, const int k, const int m);
void detStartSplitting(const double *x, const double *c, int *y, int *nums, const int n, const int m, const int k);
void blockSplitting1(const double *x, const double *c, int* const y, int* const nums, const int m, const int k, const int perProc, const int idSt);
void MPI_DetStartSplitting(const double *x, const double *c, int* const y, int* const nums, const int n, const int m, const int k);
void calcCores(const double *x, double* const c, const int* const res, const int* const nums, const int n, const int m);
void simpleCalcCores(const double* const x, double* const c, const int* const res, const int* const nums, const int m, const int perProc, const int id);
void MPI_CalcCores(const double *x, double *c, const int *res, const int *nums, const int n, const int m, const int k);
char checkSplitting(const double *x, const double *c, int *res, int *nums, const int n, const int m, const int k);
int blockSplitting2(const double *x, const double *c, int* const y, int* const nums, const int m, const int k, const int perProc, const int idSt);
char MPI_CheckSplitting(const double *x, const double *c, int* const y, int* const nums, const int n, const int m, const int k);
char constr(const int *y, const int val, const int s);
void startCoreNums(int *y, const int k, const int n);
void MPI_kmeans(const double* const X, int* const y, const int n, const int m, const int k);


#endif
