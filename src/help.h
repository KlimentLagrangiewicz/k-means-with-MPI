#ifndef HELP_H_
#define HELP_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void fscanfData(const char *fn, double *x, const int n);
void fscanfSplitting(const char *fn, int *y, const int n);
double getPrecision(int *x, int *y, const int n);
void fprintfTime(const char *fn, const double t);
void fprintfResults(const char *fn, const int *res, const int n, const int m, const int k);

#endif
