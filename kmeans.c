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

int getCluster(const double *x, const double *c, const int m, const int k, const int l) {
	int i, res = 0;
	double cur, dis = DBL_MAX;
	const double *p1 = x, *p2 = c;
	for (i = 0; i < k; i++) {
		cur = distEv(p1 + l, p2 + i * m, m);
		if (cur < dis) {
			dis = cur;
			res = i;
		}
	}
	return res;
}

static short constr(int *y, const int val, const int s) {
	int i;
	for (i = 0; i < s; i++) {
		if (y[i] == val) return 1;
	}
	return 0;
}

void startCoreNums(int *y, const int k, const int n) {
	srand((unsigned int)time(NULL));
	int i, val;
	for (i = 0; i < k; i++) {
		do {
			val = rand() % n;
		} while (constr(y, val, i));
		y[i] = val;
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
