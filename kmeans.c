#include "kmeans.h"

double distEv(const double *x, const double *c, const int m) {
	double d, r = 0;
	int i = 0;
	while (i++ < m) {
		d = *(x++) - *(c++);
		r += d * d;
	}
	return r;
}

int getCluster(const double *x, const double *c, const int m, const int k) {
	double curD, minD = DBL_MAX;
	int counter, res;
	counter = res = 0;
	while (counter < k) {
		curD = distEv(x, c, m);
		if (curD < minD) {
			minD = curD;
			res = counter;
		}
		counter++;
		c += m;
	}
	return res;
}

static short constr(const int *y, const int val, const int s) {
	int i = 0;
	while (i < s) {
		if (*(y++) == val) return 1;
		i++;
	}
	return 0;
}

void startCoreNums(int *y, const int k, const int n) {
	srand((unsigned int)time(NULL));
	int i = 0, val;
	while (i < k) {
		do {
			val = rand() % n;
		} while (constr(&y[0], val, i));
		y[i] = val;
		i++;
	}
}

void normalization(double *x, const int n, const int m) {
	int i, j;
	double mean, d, s;
	for (j = 0; j < m; j++) {
		mean = s = 0;
		for (i = j; i < n * m; i += m) {
			mean += x[i];
		}
		mean /= n;
		for (i = j; i < n * m; i += m) {
			d = x[i] - mean;
			x[i] = d;
			s += d * d;

		}
		s = sqrt(s / n);
		for (i = j; i < n * m; i += m) {
			x[i] /= s;

		}
	}
}
