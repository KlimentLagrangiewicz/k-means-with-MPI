#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"
#include "help.h"
#include "kmeans.h"

int main(int argc, char **argv) {
	if (argc < 6) {
		printf("Not enough parameters...\n");
		exit(1);
	}
	int n = atoi(argv[1]), m = atoi(argv[2]), k = atoi(argv[3]), rank, numOfProc, i, j, buf1, buf2, count, countbuf, f;
	double q, t;
	int *r = (int*)malloc(n * sizeof(int));
	int *rbuf = (int*)malloc(n * sizeof(int));
	double *x = (double*)malloc(n * m * sizeof(double));
	double *c = (double*)malloc(k * m * sizeof(double));
	double *nc = (double*)malloc(k * m * sizeof(double));
	memset(nc, 0, k * m * sizeof(double));
	int *numsBuf = (int*)malloc(k * sizeof(int));
	int *nums = (int*)malloc(k * sizeof(int));
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		fscanfData(x, n * m, argv[4]);
	}
	t = MPI_Wtime();
	if (rank == 0) {
		startCoreNums(nums, k, n);
		normalization(x, n, m);
	}
	MPI_Bcast(x, n * m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(nums, k, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProc);
	const short flag = (numOfProc > 2 * k) ? 1 : 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (flag) {
    	for (i = rank; i < k * m; i += numOfProc) {
    		nc[i] = x[nums[i / m] * m + i % m];
    	}
	} else {
		for (i = rank; i < k; i += numOfProc) {
			buf1 = nums[i] * m;
			buf2 = i * m;
			for (j = 0; j < m; j++) {
				nc[buf2 + j] = x[buf1 + j];
			}
		}
	}
	MPI_Allreduce(nc, c, k * m, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	do {
		countbuf = 0;
		memset(rbuf, 0, n * sizeof(int));
		memset(numsBuf, 0, k * sizeof(int));
		memset(nc, 0, k * m * sizeof(double));
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		for (i = rank; i < n; i += numOfProc) {
			buf2 = i * m;
			f = getCluster(x, c, m, k, buf2);
			numsBuf[f]++;
			if (f == r[i]) countbuf++;
			rbuf[i] = f;
			buf1 = f * m;
			for (j = 0; j < m; j++) {
				nc[buf1 + j] += x[buf2 + j];
			}
		}
		MPI_Allreduce(&countbuf, &count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(nc, c, k * m, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(rbuf, r, n, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(numsBuf, nums, k, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		memset(nc, 0, k * m * sizeof(double));
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		for (i = rank; i < k; i += numOfProc) {
			buf1 = nums[i];
			for (j = i * m; j < (i + 1) * m; j++) {
				nc[j] = c[j] / buf1;
			}
		}
		MPI_Allreduce(nc, c, k * m, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	} while (count != n);
	t = MPI_Wtime() - t;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		q = calcQualityOfSplitting(x, r, n, m);
		if (argc > 6) {
			int *id = (int*)malloc(n * sizeof(int));
			fscanfIdealSpliting(id, n, argv[6]);
			double a = caclAccuracy(id, r, n);
			fprintfFullResult(r, n, t, a, q, argv[5]);
			printf("Time for k-Means with MPI clustering = %lf;\nAccuracy of k-Means with MPI clustering = %lf;\nQuality of k-Means with MPI clustering = %lf;\nThe work of the program is completed!\n", t, a, q);
			free(id);
		} else {
			fprintfResult(r, n, t, q, argv[5]);
			printf("Time for k-Means with MPI clustering = %lf;\nQuality of k-Means with MPI clustering = %lf;\nThe work of the program is completed!\n", t, q);
		}
	}
	MPI_Finalize();
	free(r);
	free(rbuf);
	free(x);
	free(c);
	free(nc);
	free(numsBuf);
	free(nums);
	return 0;
}
