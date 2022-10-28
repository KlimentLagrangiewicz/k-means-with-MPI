#include "help.h"

void fscanfData(double *x, const int n, const char *fn) {
	FILE *fl;
	if ((fl = fopen(fn, "r")) == NULL) {
		printf("Error during opening data %s file...\n", fn);
		exit(1);
	}
	int i;
	for (i = 0; i < n && !feof(fl); i++) {
		if (fscanf(fl, "%lf", &x[i]) == 0) {}
	}
	fclose(fl);
}

void fprintfResult(const int *r, const int n, const double t, const double q, const char *fn) {
	FILE *fl;
	if ((fl = fopen(fn, "a")) == NULL) {
		printf("Error during opening result %s file...\n", fn);
		exit(1);
	}
	fprintf(fl, "Time for k-Means with MPI clustering = %lf;\nQuality of k-Means with MPI clustering = %lf;\n", t, q);
	int i;
	for (i = 0; i < n; i++) {
		fprintf(fl, "Object[%d]: %d;\n", i, r[i]);
	}
	fprintf(fl, "\n");
	fclose(fl);
}

void fscanfIdealSpliting(int *id, const int n, const char *fn) {
	FILE *fl;
	if ((fl = fopen(fn, "r")) == NULL) {
		printf("Error during opening spiting %s file...\n", fn);
		exit(1);
	}
	int i;
	for (i = 0; i < n && !feof(fl); i++) {
		if (fscanf(fl, "%d", &id[i]) == 0) {}
	}
	fclose(fl);
}

int getNumOfClass(const int *y, const int n) {
	int i, j, cur;
	short *v = (short*)malloc(n * sizeof(short));
	memset(v, 0, n * sizeof(short));
	for (i = 0; i < n; i++) {
		while ((v[i]) && (i < n)) i++;
		cur = y[i];
		for (j = i + 1; j < n; j++) {
			if (y[j] == cur)
				v[j] = 1;
		}
	}
	cur = 0;
	for (i = 0; i < n; i++) {
		if (v[i] == 0) cur++;
	}
	free(v);
	return cur;
}

double getCurAccuracy(const int *x, const int *y, const int *a, const int n) {
	int i, j = 0;
	for (i = 0; i < n; i++) {
		if (x[i] != a[y[i]]) j++;
	}
	return 1.0 - ((double)j / (double)n);
}

void solve(const int *x, const int *y, int *items, int size, int l, const int n, double *eps) {
    int i;
    if (l == size) {
    	double cur = getCurAccuracy(x, y, items, n);
    	if (cur > *eps) *eps = cur;
    } else {
        for (i = l; i < size; i++) {
            if (l ^ i) {
            	items[l] ^= items[i];
            	items[i] ^= items[l];
            	items[l] ^= items[i];
            	solve(x, y, items, size, l + 1, n, eps);
            	items[l] ^= items[i];
            	items[i] ^= items[l];
            	items[l] ^= items[i];
            } else {
            	solve(x, y, items, size, l + 1, n, eps);
            }
        }
    }
}

double caclAccuracy(const int *ideal, const int *r, const int n) {
	int i, j = 0, k = getNumOfClass(ideal, n);
	int *nums = (int*)malloc(k * sizeof(int));
	for (i = 0; i < k; i++) {
		nums[i] = i;
	}
	double max = getCurAccuracy(r, ideal, nums, n);
	solve(r, ideal, nums, k, j, n, &max);
	free(nums);
	return max;
}

void fprintfFullResult(const int *r, const int n, const double t, const double a, const double q, const char *fn) {
	FILE *fl;
	if ((fl = fopen(fn, "a")) == NULL) {
		printf("Error during opening result %s file...\n", fn);
		exit(1);
	}
	fprintf(fl, "Time for k-Means with MPI clustering = %lf;\nAccuracy of k-Means with MPI clustering = %lf;\nQuality of k-Means with MPI clustering = %lf;\n", t, a, q);
	int i;
	for (i = 0; i < n; i++) {
		fprintf(fl, "Object[%d]: %d;\n", i, r[i]);
	}
	fprintf(fl, "\n");
	fclose(fl);
}

static double distance(const double *x, const double *c, const int m) {
	double d, r = 0;
	int i = 0;
	while (i++ < m) {
		d = *(x++) - *(c++);
		r += d * d;
	}
	return sqrt(r);
}

double calcQualityOfSplitting(const double *x, const int *r, const int n, const int m) {
	int i, j, k1 = 0, k2 = 0, buf;
	double inside, outside;
	inside = outside = 0;
	for (i = 0; i < n; i++) {
		buf = i * m;
		for (j = i + 1; j < n; j++)  {
			if (r[i] == r[j]) {
				inside += distance(&x[buf], &x[j * m], m);
				k1++;
			} else {
				outside += distance(&x[buf], &x[j * m], m);
				k2++;
			}
		}
	}
	inside /= (k1 == 0) ? 1 : k1;
	outside /= (k2 == 0) ? 1 : k2;
	return (inside / ((outside == 0) ? inside : outside));
}
