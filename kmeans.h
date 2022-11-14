#ifndef KMEANS_H_
#define KMEANS_H_

#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <math.h>

double distEv(const double *x, const double *c, const int m);
int getCluster(const double *x, const double *c, const int m, const int k);
void startCoreNums(int *y, const int k, const int n);
void normalization(double *x, const int n, const int m);

#endif
