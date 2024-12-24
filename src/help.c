#include "help.h"


void fscanfData(const char *fn, double *x, const int n) {
	FILE *fl = fopen(fn, "r");
	if (!fl) {
		printf("Error in opening %s file...\n", fn);
		exit(1);
	}
	int i = 0;
	while (i < n && !feof(fl)) {
		if (fscanf(fl, "%lf", x + i) == 0) {}
		i++;
	}
	fclose(fl);
}

void fscanfSplitting(const char *fn, int *y, const int n) {
	FILE *fl = fopen(fn, "r");
	if (!fl) {
		printf("Can't access %s file with ideal splitting for reading...\n", fn);
		exit(1);
	}
	int i = 0;
	while (i < n && !feof(fl)) {
		if (fscanf(fl, "%d", y + i) == 0) {
			printf("Error in reading the perfect partition from %s file\n", fn);
			exit(1);
		}
		i++;
	}
	fclose(fl);
}

double getPrecision(int *x, int *y, const int n) {
	int i, j, yy = 0, ny = 0;
	for (i = 0; i < n; i++) {
		for (j = i + 1; j < n; j++) {
			if (x[i] == x[j] && y[i] == y[j]) yy++;
			if (x[i] != x[j] && y[i] == y[j]) ny++;
		}
	}
	return yy == 0 && ny == 0 ? 0.0 : (double)yy / (double)(yy + ny);
}

void fprintfTime(const char *fn, const double t) {
	FILE *fl = fopen(fn, "a");
	if (!fl) {
		printf("Error in opening %s file...\n", fn);
		exit(1);
	}
	fprintf(fl, "%lf\n", t);
	fclose(fl);
}

void fprintfResults(const char *fn, const int *res, const int n, const int m, const int k) {
	FILE *fl = fopen(fn, "a");
	if (!fl) {
		printf("Error in opening %s file...\n", fn);
		exit(1);
	}
	fprintf(fl, "Results of clustering using k-means...\nParameters: n = %d, m = %d, k = %d;\n", n, m, k);
	int i = 0;
	while (i < n) {
		fprintf(fl, "Object [%d] = %d;\n", i, res[i]);
		i++;
	}
	fclose(fl);
}
