#include "kmeans.h"

double getEvDist(const double *x1, const double *x2, const int m) {
	double d, r = 0;
	int i = 0;
	while (i++ < m) {
		d = *(x1++) - *(x2++);
		r += d * d;
	}
	return sqrt(r);
}

void autoscaling(double* const x, const int n, const int m) {
	const int s = n * m;
	double sd, Ex, Exx;
	int i, j = 0;
	while (j < m) {
		i = j;
		Ex = Exx = 0;
		while (i < s) {
			sd = x[i];
			Ex += sd;
			Exx += sd * sd;
			i += m;
		}
		Exx /= n;
		Ex /= n;
		sd = sqrt(Exx - Ex * Ex);
		i = j;
		while (i < s) {
			x[i] = (x[i] - Ex) / sd;
			i += m;
		}
		j++;
	}
}

void blockFunction1(const double* const x, double* const Ex, double* const Exx, const int m, const int perProc, const int rowID) {
	double val;
	int i, j, k = rowID * m;
	for (i = 0; i < perProc; i++) {
		for (j = 0; j < m; j++) {
			val = x[k];
			Ex[j] += val;
			Exx[j] += val * val;
			k++;
		}
	}
}

void blockFunction2(double* const x, const double* const Ex, const double* const Exx, const int m, const int perProc, const int rowID) {
	int i, j, k = rowID * m;
	for (i = 0; i < perProc; i++) {
		for (j = 0; j < m; j++) {
			x[k] = (x[k] - Ex[j]) / Exx[j];
			k++;
		}
	}
}

void MPI_Autoscaling(double* const x, const int n, const int m) {
	int pid, perProc, numOfProc;
	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProc);
	perProc = n / numOfProc;
	if (perProc == 0) {
		MPI_Comm_rank(MPI_COMM_WORLD, &pid);
		if (pid == 0) {
			autoscaling(x, n, m);
		}
		MPI_Bcast(&x[0], n * m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	} else {
		double *Ex = (double*)malloc(m * sizeof(double));
		double *Exx = (double*)malloc(m * sizeof(double));
		memset(&Ex[0], 0, m * sizeof(double));
		memset(&Exx[0], 0, m * sizeof(double));
		MPI_Comm_rank(MPI_COMM_WORLD, &pid);
		if (pid == 0) {
			blockFunction1(x, Ex, Exx, m, perProc, 0);
			if (n > perProc * numOfProc) {
				blockFunction1(x, Ex, Exx, m, n - perProc * numOfProc, perProc * numOfProc);
			}
			int i, j;
			double *buffEx = (double*)malloc(m * sizeof(double));
			double *buffExx = (double*)malloc(m * sizeof(double));
			for (i = 1; i < numOfProc; i++) {
				MPI_Recv(&buffEx[0], m, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &status);
				MPI_Recv(&buffExx[0], m, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &status);
				for (j = 0; j < m; j++) {
					Ex[j] += buffEx[j];
					Exx[j] += buffExx[j];
				}
			}
			for (i = 0; i < m; i++) {
				Ex[i] /= n;
				Exx[i] /= n;
				Exx[i] = sqrt(Exx[i] - Ex[i] * Ex[i]);
			}
			for (i = 1; i < numOfProc; i++) {
				MPI_Send(&Ex[0], m, MPI_DOUBLE, i, 4, MPI_COMM_WORLD);
				MPI_Send(&Exx[0], m, MPI_DOUBLE, i, 5, MPI_COMM_WORLD);
			}
			blockFunction2(x, Ex, Exx, m, perProc, 0);
			if (n > perProc * numOfProc) {
				blockFunction2(x, Ex, Exx, m, n - perProc * numOfProc, perProc * numOfProc);
			}
			for (i = 1; i < numOfProc; i++) {
				MPI_Recv(&x[i * perProc * m], perProc * m, MPI_DOUBLE, i, 6, MPI_COMM_WORLD, &status);
			}
			for (i = 1; i < numOfProc; i++) {
				MPI_Send(&x[0], n * m, MPI_DOUBLE, i, 7, MPI_COMM_WORLD);
			}
			free(buffEx);
			free(buffExx);
		} else {
			blockFunction1(x, Ex, Exx, m, perProc, pid * perProc);
			MPI_Send(&Ex[0], m, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
			MPI_Send(&Exx[0], m, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
			MPI_Recv(&Ex[0], m, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD, &status);
			MPI_Recv(&Exx[0], m, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD, &status);
			blockFunction2(x, Ex, Exx, m, perProc, pid * perProc);
			MPI_Send(&x[pid * perProc * m], perProc * m, MPI_DOUBLE, 0, 6, MPI_COMM_WORLD);
			MPI_Recv(&x[0], n * m, MPI_DOUBLE, 0, 7, MPI_COMM_WORLD, &status);
		}
		free(Ex);
		free(Exx);
	}
}

int getCluster(const double *x, const double *c, const int m, const int k) {
	double curD, minD = DBL_MAX;
	int counter, res;
	counter = res = 0;
	while (counter < k) {
		curD = getEvDist(x, c, m);
		if (curD < minD) {
			minD = curD;
			res = counter;
		}
		counter++;
		c += m;
	}
	return res;
}


void detCores(const double *x, double *c, const int *sn, const int k, const int m) {
	int i;
	for (i = 0; i < k; i++) {
		memcpy(&c[i * m], &x[sn[i] * m], m * sizeof(double));
	}
}

void detStartSplitting(const double *x, const double *c, int *y, int *nums, const int n, const int m, const int k) {
	int i = 0, j = 0, cur;
	while (i < n) {
		cur = getCluster(&x[j], &c[0], m, k);
		y[i] = cur;
		nums[cur]++;
		j += m;
		i++;
	}
}

void blockSplitting1(const double *x, const double *c, int* const y, int* const nums, const int m, const int k, const int perProc, const int idSt) {
	int i, cur;
	for (i = idSt; i < idSt + perProc; i++) {
		cur = getCluster(&x[i * m], &c[0], m, k);
		y[i] = cur;
		nums[cur]++;
	}
}

void MPI_DetStartSplitting(const double *x, const double *c, int* const y, int* const nums, const int n, const int m, const int k) {
	int pid, perProc, numOfProc;
	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProc);
	perProc = n / numOfProc;
	memset(&nums[0], 0, k * sizeof(int));
	if (perProc == 0) {
		MPI_Comm_rank(MPI_COMM_WORLD, &pid);
		if (pid == 0) {
			detStartSplitting(x, c, y, nums, n, m, k);
		}
		MPI_Bcast(&y[0], n, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&nums[0], k, MPI_INT, 0, MPI_COMM_WORLD);
	} else {
		MPI_Comm_rank(MPI_COMM_WORLD, &pid);
		if (pid == 0) {
			blockSplitting1(x, c, y, nums, m, k, perProc, 0);
			if (n > numOfProc * perProc) {
				blockSplitting1(x, c, y, nums, m, k, n - numOfProc * perProc, numOfProc * perProc);
			}
			int *nums_buff = (int*)malloc(k * sizeof(int));
			int i, j;
			for (i = 1; i < numOfProc; i++) {
				MPI_Recv(&y[i * perProc], perProc, MPI_INT, i, 10, MPI_COMM_WORLD, &status);
				MPI_Recv(&nums_buff[0], k, MPI_INT, i, 11, MPI_COMM_WORLD, &status);
				for (j = 0; j < k; j++) {
					nums[j] += nums_buff[j];
				}
			}
			for (i = 1; i < numOfProc; i++) {
				MPI_Send(&y[0], n, MPI_INT, i, 12, MPI_COMM_WORLD);
				MPI_Send(&nums[0], k, MPI_INT, i, 13, MPI_COMM_WORLD);
			}
			free(nums_buff);
		} else {
			blockSplitting1(x, c, y, nums, m, k, perProc, pid * perProc);
			MPI_Send(&y[pid * perProc], perProc, MPI_INT, 0, 10, MPI_COMM_WORLD);
			MPI_Send(&nums[0], k, MPI_INT, 0, 11, MPI_COMM_WORLD);

			MPI_Recv(&y[0], n, MPI_INT, 0, 12, MPI_COMM_WORLD, &status);
			MPI_Recv(&nums[0], k, MPI_INT, 0, 13, MPI_COMM_WORLD, &status);
		}
	}
}

void calcCores(const double *x, double* const c, const int* const res, const int* const nums, const int n, const int m) {
	int i, j, buf1, buf2, buf3;
	for (i = 0; i < n; i++) {
		buf1 = nums[res[i]];
		buf2 = res[i] * m;
		buf3 = i * m;
		for (j = 0; j < m; j++) {
			c[buf2 + j] += x[buf3 + j] / buf1;
		}
	}
}

void simpleCalcCores(const double* const x, double* const c, const int* const res, const int* const nums, const int m, const int perProc, const int id) {
	int i, j;
	for (i = id; i < id + perProc; i++) {
		const int buf1 = nums[res[i]], buf2 = res[i] * m, buf3 = i * m;
		for (j = 0; j < m; j++) {
			c[buf2 + j] += x[buf3 + j] / buf1;
		}
	}
}

void MPI_CalcCores(const double *x, double *c, const int *res, const int *nums, const int n, const int m, const int k) {
	memset(&c[0], 0, k * m * sizeof(double));
	int pid, perProc, numOfProc;
	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProc);
	perProc = n / numOfProc;
	if (perProc == 0) {
		MPI_Comm_rank(MPI_COMM_WORLD, &pid);
		if (pid == 0) {
			calcCores(x, c, res, nums, n, m);
		}
		MPI_Bcast(&c[0], k * m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	} else {
		MPI_Comm_rank(MPI_COMM_WORLD, &pid);
		if (pid == 0) {
			simpleCalcCores(x, c, res, nums, m, perProc, 0);
			if (n > numOfProc * perProc) {
				simpleCalcCores(x, c, res, nums, m, n - numOfProc * perProc, numOfProc * perProc);
			}
			double *bufCores = (double*)malloc(k * m * sizeof(double));
			int i, j;
			for (i = 1; i < numOfProc; i++) {
				MPI_Recv(&bufCores[0], k * m, MPI_DOUBLE, i, 15, MPI_COMM_WORLD, &status);
				for (j = 0; j < k * m; j++) {
					c[j] += bufCores[j];
				}
			}
			for (i = 1; i < numOfProc; i++) {
				MPI_Send(&c[0], k * m, MPI_DOUBLE, i, 16, MPI_COMM_WORLD);
			}
			free(bufCores);
		} else {
			simpleCalcCores(x, c, res, nums, m, perProc, pid * perProc);
			MPI_Send(&c[0], k * m, MPI_DOUBLE, 0, 15, MPI_COMM_WORLD);

			MPI_Recv(&c[0], k * m, MPI_DOUBLE, 0, 16, MPI_COMM_WORLD, &status);
		}
	}
}

char checkSplitting(const double *x, const double *c, int *res, int *nums, const int n, const int m, const int k) {
	int i = 0, count = 0, j = 0, f;
	while (i < n) {
		f = getCluster(&x[j], &c[0], m, k);
		if (f == res[i]) count++;
		res[i] = f;
		nums[f]++;
		j += m;
		i++;
	}
	return (n == count) ? 0 : 1;
}

int blockSplitting2(const double *x, const double *c, int* const y, int* const nums, const int m, const int k, const int perProc, const int idSt) {
	int i, cur, res = 0;
	for (i = idSt; i < idSt + perProc; i++) {
		cur = getCluster(&x[i * m], &c[0], m, k);
		if (cur == y[i]) res++;
		y[i] = cur;
		nums[cur]++;
	}
	return res;
}

char MPI_CheckSplitting(const double *x, const double *c, int* const y, int* const nums, const int n, const int m, const int k) {
	int pid, perProc, numOfProc;
	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProc);
	perProc = n / numOfProc;
	char result = 1;
	memset(&nums[0], 0, k * sizeof(int));
	if (perProc == 0) {
		MPI_Comm_rank(MPI_COMM_WORLD, &pid);
		if (pid == 0) {
			result = checkSplitting(x, c, y, nums, n, m, k);
		}
		MPI_Bcast(&result, 1, MPI_CHAR, 0, MPI_COMM_WORLD);
		MPI_Bcast(&nums[0], k, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&y[0], n, MPI_INT, 0, MPI_COMM_WORLD);
	} else {
		MPI_Comm_rank(MPI_COMM_WORLD, &pid);
		if (pid == 0) {
			int count = blockSplitting2(x, c, y, nums, m, k, perProc, 0), bufCount;
			if (n > numOfProc * perProc) {
				count += blockSplitting2(x, c, y, nums, m, k, n - numOfProc * perProc, numOfProc * perProc);
			}
			int *numsBuff = (int*)malloc(k * sizeof(int));
			int i, j;
			for (i = 1; i < numOfProc; i++) {
				MPI_Recv(&bufCount, 1, MPI_INT, i, 20, MPI_COMM_WORLD, &status);
				count += bufCount;
				MPI_Recv(&y[i * perProc], perProc, MPI_INT, i, 21, MPI_COMM_WORLD, &status);
				MPI_Recv(&numsBuff[0], k, MPI_INT, i, 22, MPI_COMM_WORLD, &status);
				for (j = 0; j < k; j++) {
					nums[j] += numsBuff[j];
				}
			}
			free(numsBuff);
			result = (count == n) ? 0 : 1;
			for (i = 1; i < numOfProc; i++) {
				MPI_Send(&result, 1, MPI_CHAR, i, 23, MPI_COMM_WORLD);
				MPI_Send(&nums[0], k, MPI_INT, i, 24, MPI_COMM_WORLD);
				MPI_Send(&y[0], n, MPI_INT, i, 25, MPI_COMM_WORLD);
			}
		} else {
			int count = blockSplitting2(x, c, y, nums, m, k, perProc, pid * perProc);
			MPI_Send(&count, 1, MPI_INT, 0, 20, MPI_COMM_WORLD);
			MPI_Send(&y[pid * perProc], perProc, MPI_INT, 0, 21, MPI_COMM_WORLD);
			MPI_Send(&nums[0], k, MPI_INT, 0, 22, MPI_COMM_WORLD);

			MPI_Recv(&result, 1, MPI_CHAR, 0, 23, MPI_COMM_WORLD, &status);
			MPI_Recv(&nums[0], k, MPI_INT, 0, 24, MPI_COMM_WORLD, &status);
			MPI_Recv(&y[0], n, MPI_INT, 0, 25, MPI_COMM_WORLD, &status);
		}
	}
	return result;
}

char constr(const int *y, const int val, const int s) {
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

void MPI_kmeans(const double* const X, int* const y, const int n, const int m, const int k) {
	double *x = (double*)malloc(n * m * sizeof(double));
	int pid;
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	if (pid == 0) {
		memcpy(&x[0], &X[0], n * m * sizeof(double));
	}
	MPI_Bcast(&x[0], n * m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Autoscaling(x, n, m);
	int *nums = (int*)malloc(k * sizeof(double));
	double *c = (double*)malloc(k * m * sizeof(double));
	if (pid == 0) {
		startCoreNums(nums, k, n);
		detCores(x, c, nums, k, m);
	}
	MPI_Bcast(&c[0], k * m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_DetStartSplitting(x, c, y, nums, n, m, k);
	char flag = 1;
	do {
		MPI_CalcCores(x, c, y, nums, n, m, k);
		flag = MPI_CheckSplitting(x, c, y, nums, n, m, k);
		MPI_Barrier(MPI_COMM_WORLD);
	} while (flag);
	free(c);
	free(x);
	free(nums);
}
