#include "kmeans.h"


double getEvDist(const double *x1, const double *x2, int m) {
	double d, r = 0;
	while (m--) {
		d = *(x1++) - *(x2++);
		r += d * d;
	}
	return sqrt(r);
}

void autoscaling(double* const x, const int n, const int m) {
	const int s = n * m;
	int j;
	for (j = 0; j < m; j++) {
		double sd, Ex = 0.0, Exx = 0.0, *ptr;
		for (ptr = x + j; ptr < x + s; ptr += m) {
			sd = *ptr;
			Ex += sd;
			Exx += sd * sd;
		}
		Exx /= n;
		Ex /= n;
		sd = sqrt(Exx - Ex * Ex);
		for (ptr = x + j; ptr < x + s; ptr += m) {
			*ptr = (*ptr - Ex) / sd;
		}
	}
}

void blockFunction1(const double* const x, double* const Ex, double* const Exx, const int m, const int perProc) {
	int i, j;
	for (i = 0; i < perProc * m; i += m) { 
		for (j = 0; j < m; j++) {
			Ex[j] += x[i + j];
			Exx[j] += x[i + j] * x[i + j];
		}
	}
}

void blockFunction2(double* const x, const double* const Ex, const double* const Exx, const int m, const int perProc) {
	int i, j;
	for (i = 0; i < perProc * m; i += m)
		for (j = 0; j < m; j++) 
			x[i + j] = (x[i + j] - Ex[j]) / Exx[j];
}

void MPI_Autoscaling(double* const x, const int n, const int m) {
	int numOfProc;
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProc);
	int perProc = n / numOfProc;
	if (perProc == 0) {
		int pid;
		MPI_Comm_rank(MPI_COMM_WORLD, &pid);
		if (pid == 0) autoscaling(x, n, m);
	} else {
		double *localX = (double*)malloc(perProc * m * sizeof(double));
		MPI_Scatter(x, perProc * m, MPI_DOUBLE, localX, perProc * m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		double *Ex = (double*)calloc(m, sizeof(double));
		double *Exx = (double*)calloc(m, sizeof(double));
		int pid;
		MPI_Comm_rank(MPI_COMM_WORLD, &pid);
		if (pid == 0 && n > perProc * numOfProc) blockFunction1(x + perProc * numOfProc * m, Ex, Exx, m, n - perProc * numOfProc);
		blockFunction1(localX, Ex, Exx, m, perProc);
		if (pid == 0) MPI_Reduce(MPI_IN_PLACE, Ex, m, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		else MPI_Reduce(Ex, Ex, m, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if (pid == 0) MPI_Reduce(MPI_IN_PLACE, Exx, m, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		else MPI_Reduce(Exx, Exx, m, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if (pid == 0) {
			int i;
			for (i = 0; i < m; i++) {
				Ex[i] /= n;
				Exx[i] /= n;
				Exx[i] = sqrt(Exx[i] - Ex[i] * Ex[i]);
			}
		}
		MPI_Bcast(Ex, m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(Exx, m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if (pid == 0 && n > perProc * numOfProc) blockFunction2(x + perProc * numOfProc * m, Ex, Exx, m, n - perProc * numOfProc);
		blockFunction2(localX, Ex, Exx, m, perProc);
		MPI_Gather(localX, perProc * m, MPI_DOUBLE, x, perProc * m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		free(Ex);
		free(Exx);
		free(localX);
	}
}

int getCluster(const double* const x, const double* const c, const int m, int k) {
	int res = --k;
	double minD = getEvDist(x, c + k * m, m);
	while (k--) {
		const double curD = getEvDist(x, c + k * m, m);
		if (curD < minD) {
			minD = curD;
			res = k;
		}
	}
	return res;
}

char constr(const int *y, const int val, int s) {
	while (s--) {
		if (*(y++) == val) return 1;
	}
	return 0;
}

void detCores(const double* const x, double* const c, const int n, const int m, const int k) {
	int *nums = (int*)malloc(k * sizeof(int));
	srand((unsigned int)clock());
	int i;
	for (i = 0; i < k; i++) {
		int val = rand() % n;
		while (constr(nums, val, i)) {
			val = rand() % n;
		}
		nums[i] = val;
		memcpy(c + i * m, x + val * m, m * sizeof(double));
	}
	free(nums);
}

void detStartPartition(const double* const x, const double* const c, int* const y, int* const nums, int n, const int m, const int k) {
	while (n--) {
		const int l = getCluster(x + n * m, c, m, k);
		y[n] = l;
		nums[l]++;
	}
}

void MPI_Detstartpartition(const double* const x, const double* const c, int* const y, int* const nums, const int n, const int m, const int k) {
	int numOfProc;
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProc);
	int perProc = n / numOfProc;
	if (perProc == 0) {
		int pid;
		MPI_Comm_rank(MPI_COMM_WORLD, &pid);
		if (pid == 0) {
			memset(nums, 0, k * sizeof(int));
			detStartPartition(x, c, y, nums, n, m, k);
		}
		MPI_Bcast(nums, k, MPI_INT, 0, MPI_COMM_WORLD);
	} else {
		double *localX = (double*)malloc(perProc * m * sizeof(double));
		MPI_Scatter(x, perProc * m, MPI_DOUBLE, localX, perProc * m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		int *localY = (int*)malloc(perProc * sizeof(int));
		memset(nums, 0, k * sizeof(int));
		detStartPartition(localX, c, localY, nums, perProc, m, k);
		int pid;
		MPI_Comm_rank(MPI_COMM_WORLD, &pid);
		if (pid == 0 && n > perProc * numOfProc) detStartPartition(x + perProc * numOfProc * m, c, y + perProc * numOfProc, nums, n - perProc * numOfProc, m, k);
		MPI_Allreduce(MPI_IN_PLACE, nums, k, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Gather(localY, perProc, MPI_INT, y, perProc, MPI_INT, 0, MPI_COMM_WORLD);
		free(localX);
		free(localY);
	}
}

void sumCores(const double* const x, double* const c, const int* const y, const int n, const int m) {
	int i, j;
	for (i = 0; i < n; i++) {
		const int f = y[i] * m, l = i * m;
		for (j = 0; j < m; j++) {
			c[f + j] += x[l + j];
		}
	}
}

void MPI_Sumcores(const double* const x, double* const c, const int* const y, const int n, const int m, const int k) {
	int numOfProc;
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProc);
	int perProc = n / numOfProc;
	if (perProc == 0) {
		int pid;
		MPI_Comm_rank(MPI_COMM_WORLD, &pid);
		memset(c, 0, k * m * sizeof(double));
		if (pid == 0) sumCores(x, c, y, n, m);
		MPI_Bcast(c, k * m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	} else {
		memset(c, 0, k * m * sizeof(double));
		double *localX = (double*)malloc(perProc * m * sizeof(double));
		MPI_Scatter(x, perProc * m, MPI_DOUBLE, localX, perProc * m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		int *localY = (int*)malloc(perProc * sizeof(int));
		MPI_Scatter(y, perProc, MPI_INT, localY, perProc, MPI_INT, 0, MPI_COMM_WORLD);
		int pid;
		MPI_Comm_rank(MPI_COMM_WORLD, &pid);
		if (pid == 0 && n > perProc * numOfProc) sumCores(x + perProc * numOfProc * m, c, y + perProc * numOfProc, n - perProc * numOfProc, m);
		sumCores(localX, c, localY, perProc, m);
		MPI_Allreduce(MPI_IN_PLACE, c, k * m, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		free(localX);
		free(localY);
	}
}

void coresDiv(double* const c, const int* const nums, const int m, const int k) {
	int i, j;
	for (i = 0; i < k; i++) {
		const int l = i * m, f = nums[i];
		for (j = l; j < l + m; j++)
			c[j] /= f;
	}
}

void MPI_Coresdiv(double* const c, const int* const nums, const int m, const int k) {
	int numOfProc;
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProc);
	int perProc = k / numOfProc;
	if (perProc == 0) {
		int pid;
		MPI_Comm_rank(MPI_COMM_WORLD, &pid);
		if (pid == 0) coresDiv(c, nums, m, k);
		MPI_Bcast(c, k * m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	} else {
		int pid;
		MPI_Comm_rank(MPI_COMM_WORLD, &pid);
		if (pid == 0 && k > numOfProc * perProc) coresDiv(c + perProc * numOfProc * m, nums + perProc * numOfProc, m, k - perProc * numOfProc);
		coresDiv(c + pid * perProc * m, nums + perProc * pid, m, perProc);
		MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, c, perProc * m, MPI_DOUBLE, MPI_COMM_WORLD);
		if (k > numOfProc * perProc) MPI_Bcast(c + perProc * numOfProc * m, k * m - perProc * numOfProc * m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
}

char checkPartition(const double* const x, const double* const c, int* const y, int* const nums, const int n, const int m, const int k) {
	char flag = 0;
	int i;
	for (i = 0; i < n; i++) {
		const int f = getCluster(x + i * m, c, m, k);
		if (y[i] != f) flag = 1;
		y[i] = f;
		nums[f]++;
	}
	return flag;
}

char MPI_Checkpartition(const double* const x, const double* const c, int* const y, int* const nums, const int n, const int m, const int k) {
	int numOfProc;
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProc);
	int perProc = n / numOfProc;
	if (perProc == 0) {
		int pid;
		MPI_Comm_rank(MPI_COMM_WORLD, &pid);
		char flag = 0;
		if (pid == 0) {
			memset(nums, 0, k * sizeof(int));
			flag = checkPartition(x, c, y, nums, n, m, k);
		}
		MPI_Bcast(nums, k, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&flag, 1, MPI_CHAR, 0, MPI_COMM_WORLD);
		return flag;
	}
	double *localX = (double*)malloc(perProc * m * sizeof(double));
	MPI_Scatter(x, perProc * m, MPI_DOUBLE, localX, perProc * m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	int *localY = (int*)malloc(perProc * sizeof(int));
	MPI_Scatter(y, perProc, MPI_INT, localY, perProc, MPI_INT, 0, MPI_COMM_WORLD);
	memset(nums, 0, k * sizeof(int));
	char flag = checkPartition(localX, c, localY, nums, perProc, m, k);
	int pid;
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);	
	if (pid == 0 && n > numOfProc * perProc) flag |= checkPartition(x + numOfProc * perProc * m, c, y + numOfProc * perProc, nums, n - numOfProc * perProc, m, k);
	MPI_Gather(localY, perProc, MPI_INT, y, perProc, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, nums, k, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &flag, 1, MPI_CHAR, MPI_BOR, MPI_COMM_WORLD);
	free(localX);
	free(localY);
	return flag;
}

char MPI_Cyclicrecalc(const double* const x, double* const c, int* const y, int* const nums, const int n, const int m, const int k) {
	MPI_Sumcores(x, c, y, n, m, k);
	MPI_Coresdiv(c, nums, m, k);
	return MPI_Checkpartition(x, c, y, nums, n, m, k);
}

void MPI_Kmeans(const double* const X, int* const y, const int n, const int m, const int k) {
	double *x = NULL;
	int pid;
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	if (pid == 0) {
		x = (double*)malloc(n * m * sizeof(double));
		memcpy(x, X, n * m * sizeof(double));
	}
	MPI_Autoscaling(x, n, m);	
	double *c = (double*)malloc(k * m * sizeof(double));	
	if (pid == 0) detCores(x, c, n, m, k);
	MPI_Bcast(c, k * m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	int *nums = (int*)malloc(k * sizeof(double));
	MPI_Detstartpartition(x, c, y, nums, n, m, k);
	while (MPI_Cyclicrecalc(x, c, y, nums, n, m, k));
	free(c);
	free(nums);
	if (pid == 0) free(x);
}