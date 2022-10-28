#ifndef HELP_H_
#define HELP_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

void fscanfData(double *x, const int n, const char *fn);
void fprintfResult(const int *r, const int n, const double t, const double q, const char *fn);
void fscanfIdealSpliting(int *id, const int n, const char *fn);
int getNumOfClass(const int *y, const int n);
double getCurAccuracy(const int *x, const int *y, const int *a, const int n);
void solve(const int *x, const int *y, int *items, int size, int l, const int n, double *eps);
double caclAccuracy(const int *ideal, const int *r, const int n);
void fprintfFullResult(const int *r, const int n, const double t, const double a, const double q, const char *fn);
double calcQualityOfSplitting(const double *x, const int *r, const int n, const int m);

#endif
