#include <stdio.h>
#include <string.h>
#include "mpi.h"

#include "help.h"
#include "kmeans.h"

int main(int argc, char* argv[]) {
	if (argc < 6) {
		puts("Not enough parameters...");
		exit(1);
	}
	MPI_Init(&argc, &argv);
	const int n = atoi(argv[2]), m = atoi(argv[3]), k = atoi(argv[4]);
	if (n < 0 || m < 0 || k < 0 || k > n) {
		puts("Value of parameters is incorrect...");
		exit(1);
	}
	double *x = NULL;
	int *y  = NULL;
	int pid;
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	if (pid == 0) {
		x = (double*)malloc(n * m * sizeof(double));
		y = (int*)malloc(n * sizeof(int));
		fscanfData(argv[1], x, n * m);
	}	
	double t1 = MPI_Wtime();
	MPI_Kmeans(x, y, n, m, k);
	t1 = MPI_Wtime() - t1;
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);	
	if (pid == 0) {
		printf("Time for k-means clustering with usage MPI technology: %lf", t1);
		fprintfResults(argv[5], y, n, m, k);
		if (argc > 6) {			
			int *ideal = (int*)malloc(n * sizeof(int));
			fscanfSplitting(argv[6], ideal, n);
			const double p = getPrecision(ideal, y, n);
			free(ideal);
			printf("\nPrecision of k-means clustering with usage MPI technology: %lf;\n", p);
		}
		free(x);
		free(y);
	}
	MPI_Finalize();
	return 0;
}
