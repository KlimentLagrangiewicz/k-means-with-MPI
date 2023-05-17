#include <stdio.h>
#include <string.h>
#include "mpi.h"

#include "help.h"
#include "kmeans.h"

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	const int n = atoi(argv[1]), m = atoi(argv[2]), k = atoi(argv[3]);
	if ((n < 0) || (m < 0) || (k < 0) || (k > n)) {
		puts("Value of parameters is incorrect...");
			exit(1);
	}
	double *x = (double*)malloc(n * m * sizeof(double));
	int *y1  = (int*)malloc(n * sizeof(int));
	int pid;
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	if (pid == 0) {
		fscanfData(argv[4], x, n * m);
	}
	MPI_Bcast(&x[0], n * m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	double t1 = MPI_Wtime();
	MPI_kmeans(x, y1, n, m, k);
	t1 = MPI_Wtime() - t1;
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	if (pid == 0) {
			printf("Time for k-means clustering with usage MPI technology: %.6lf", t1);
			fprintfResults(argv[5], y1, n, m, k);
			if (argc > 6) {
				double a;
				int *ideal = (int*)malloc(n * sizeof(int));
				fscanfSplitting(argv[6], ideal, n);
				a = getAccuracy(ideal, y1, n);
				free(ideal);
				printf("Accuracy of clustering with usage MPI technology = %.5lf;", a);
			}
	}
	free(x);
	free(y1);
	MPI_Finalize();
	return 0;
}
