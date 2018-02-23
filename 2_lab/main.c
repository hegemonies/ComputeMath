#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define num 3
#define E 0.001

double *matrxMultVect(double *m, double *v, int n)
{
	if (!m || !v) {
		return NULL;
	}

	double *tmp_vector = calloc(n, sizeof(double));

	if (!tmp_vector) {
		return NULL;
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			tmp_vector[i] += m[i * n + j] * v[j];
		}
	}

	return tmp_vector;
}

double *vectSublVect(double *v1, double *v2, int n)
{
	if (!v1 || !v2) {
		return NULL;
	}

	double *tmp_vector = calloc(n, sizeof(double));

	for (int i = 0; i < n; i++) {
		tmp_vector[i] = v1[i] - v2[i];
	}

	return tmp_vector;
}

double getMaxFromMatrix(double *m, int n)
{
	if (!m) {
		return -1;
	}

	double max = 0;
	double sum;

	for (int i = 0; i < n; i++) {
		sum = 0;

		for (int j = 0; j < n; j++) {
			sum += fabs(m[i * n + j]);
		}

		if (sum > max) {
			max = sum;
		}
	}

	return max;
}

double getMaxFromVector(double *v, int n)
{
	if (!v) {
		return -1;
	}

	double max = 0;

	for (int i = 0; i < n; i++) {
		if (v[i] > max) {
			max = v[i];
		}
	}

	return max;
}

int getN(double maxC, double maxB, double e)
{
	return ( (log((e * (1 - maxC)) / maxB)) / (log(maxC)) ) + 1;
}

int getNumberIteration(double *m, double *v, int n, double e)
{
	if (!m || !v) {
		return -1;
	}

	double maxC = getMaxFromMatrix(m, n);
	double maxB = getMaxFromVector(v, n);

	return getN(maxC, maxB, e);
}

void redMatxToConvForm(double *m, double *v, int n)
{
	if (!m || !v) {
		return;
	}

	for (int i = 0; i < n; i++) {
		double div = m[i * n + i];
		for (int j = 0; j < n; j++) {
			m[i * n + j] /= div;
		}
		v[i] /= div;
	}
}

void zeroMainDiagonal(double *m, int n)
{
	if (!m) {
		return;
	}

	for (int i = 0; i < n; i++) {
		m[i * n + i] = 0;
	}
}

void zeroingVectors(double *v, int n)
{
	if (!v) {
		return;
	}

	for (int i = 0; i < n; i++) {
		v[i] = 0;
	}
}

void swapVectors(double *v1, double *v2, int n)
{
	if (!v1 || !v2) {
		return;
	}

	double tmp;

	for (int i = 0; i < n; i++) {
		tmp = v1[i];
		v1[i] = v2[i];
		v2[i] = tmp;
	}
}

double *methodSimpleIteration(double *A, double *B, int n)
{
	redMatxToConvForm(A, B, num); // привидение матрица к удобному виду

	zeroMainDiagonal(A, num); // зануление элементов главной диагонали

	int N = getNumberIteration(A, B, num, E); // подсчет кол-ва итераций

	printf("%d\n", N);

	double **X = calloc(2, sizeof(double));
	for (int i = 0; i < 2; i++) {
		X[i] = calloc(num, sizeof(double));
	}

	double *Cx;

	for (int i = 0; i < N; i++) {
		Cx = matrxMultVect(A, X[0], num);
		X[1] = vectSublVect(B, Cx, num);
		swapVectors(X[0], X[1], num);
		zeroingVectors(X[1], num);
	}

	return X[0];
}

double *methodSeidels(double *A, double *B, int n)
{
	redMatxToConvForm(A, B, num); // привидение матрица к удобному виду

	zeroMainDiagonal(A, num); // зануление элементов главной диагонали

	int N = getNumberIteration(A, B, num, E); // подсчет кол-ва итераций

	double **X = calloc(N, sizeof(double)); // k
	for (int i = 0; i < N; i++) {
		X[i] = calloc(num, sizeof(double));
	}

	double Cx;

	for (int k = 0; k < N - 1; k++) {
		for (int i = 0; i < n; i++) {
			for (int j1 = 0; j1 < i; j1++) {
				Cx += X[k + 1][j1] * A[i * n + j1];
			}
			for (int j = i; j < n; j++) {
				Cx += X[k][j] * A[i * n + j];
			}

			X[k + 1][i] = B[i] - Cx;
			Cx = 0;
		}
	}

	return X[N - 1];
}

int main()
{
	double A[][num] = {
		{5, -1, -1},
		{-1, -3, 0}, 
		{1, 1, 4}
	};

	double B[num] = {
		3,
		-7,
		3
	};

	// double *x = methodSimpleIteration(A, B, num);
	

	// for (int i = 0; i < num; i++) {
	// 	printf("%0.3f ", x[i]);
	// }

	// printf("\n");

	double *x1 = methodSeidels(A, B, num);
	

	for (int i = 0; i < num; i++) {
		printf("%0.3f ", x1[i]);
	}

	printf("\n");

	return 0;
}