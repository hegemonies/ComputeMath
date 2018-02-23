#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#define num 3
#define E 0.001

double *matrxMultVect(double **m, double *v, int n)
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
			tmp_vector[i] += m[i][j] * v[j];
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

double getMaxFromMatrix(double **m, int n)
{
	if (!m) {
		return -1;
	}

	double max = 0;
	double sum;

	for (int i = 0; i < n; i++) {
		sum = 0;

		for (int j = 0; j < n; j++) {
			sum += fabs(m[i][j]);
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

double getN(double maxC, double maxB, double e)
{
	return (( (log((e * (1 - maxC)) / maxB)) / (log(maxC)) ) + 1);
}

int getNumberIteration(double **m, double *v, int n, double e)
{
	if (!m || !v) {
		return -1;
	}

	double maxC = getMaxFromMatrix(m, n);

	if (maxC >= 1) {
		printf("Не сходится\n");
		exit(1);
	}

	double maxB = getMaxFromVector(v, n);

	return getN(maxC, maxB, e);
}

void redMatxToConvForm(double **m, double *v, int n)
{
	if (!m || !v) {
		return;
	}

	for (int i = 0; i < n; i++) {
		double div = m[i][i];
		if (div == 0) {
			continue;
		}
		for (int j = 0; j < n; j++) {
			m[i][j] /= div;
		}
		v[i] /= div;
	}
}

void zeroMainDiagonal(double **m, int n)
{
	if (!m) {
		return;
	}

	for (int i = 0; i < n; i++) {
		m[i][i] = 0;
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

double *methodSimpleIteration(double **A, double *B, int n)
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

double *methodSeidels(double **A, double *B, int n)
{
	redMatxToConvForm(A, B, num); // привидение матрица к удобному виду

	zeroMainDiagonal(A, num); // зануление элементов главной диагонали

	int N = getNumberIteration(A, B, num, E); // подсчет кол-ва итераций

	printf("%d\n", N);

	double **X = calloc(N, sizeof(double)); // k
	for (int i = 0; i < N; i++) {
		X[i] = calloc(num, sizeof(double));
	}

	double Cx;

	for (int k = 0; k < N - 1; k++) {
		for (int i = 0; i < n; i++) {
			for (int j1 = 0; j1 < i; j1++) {
				Cx += X[k + 1][j1] * A[i][j1];
			}
			for (int j = i; j < n; j++) {
				Cx += X[k][j] * A[i][j];
			}

			X[k + 1][i] = B[i] - Cx;
			Cx = 0;
		}
	}

	return X[N - 1];
}
/*--------------------------------------------------------*/
double wtime()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec * 1E-6;
}

int getSheetFromFile(char *file_name, double **m, double *v, int *h, int *w)
{
	FILE *in = fopen(file_name, "r");

	if (!in) {
		return -1;
	}

	char *str = NULL;
	size_t len = 0;
	*h = 0;
	*w = 0;
	int count;

	while (getline(&str, &len, in) != -1) {
		count = 0;
		for (int i = 0; str[i] != 0; i++) {
			if (str[i] == ' ') {
				count++;
			}
		}
		count++;
		if (count > *w) {
			*w = count;
		}
		(*h)++;
	}

	fseek(in, 0, SEEK_SET);

	m = realloc(m, *h);
	v = realloc(v, *h);

	for (int i = 0; i < *h; i++) {
		m[i] = calloc(*w, sizeof(double));
		for (int j = 0; j < *w - 1; j++) {
			fscanf(in, "%lf", &m[i][j]);
		}
		fscanf(in, "%lf", &v[i]);
	}

	fclose(in);

	return 0;
}


void printMatx(double **m, int n)
{
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%.3f ", m[i][j]);
		}
		printf("\n");
	}
}

void printVector(double *v, int n)
{
	for (int i = 0; i < n; i++) {
		printf("%.3f ", v[i]);
	}

	printf("\n");
}
/*----------------------------------------------------------*/

int main()
{
	// double A[][num] = {
	// 	{5, -1, -1},
	// 	{-1, -3, 0}, 
	// 	{1, 1, 4}
	// };

	// double B[num] = {
	// 	3,
	// 	-7,
	// 	3
	// };

	double **A = malloc(0);
	double *B = malloc(0);
	int height = 0;
	int width = 0;
	getSheetFromFile("in.txt", A, B, &height, &width);

	printMatx(A, num);
	printVector(B, num);

	// srand(time(0));

	// double **A = calloc(num, sizeof(double));
	// for (int i = 0; i < num; i++) {
	// 	A[i] = calloc(num, sizeof(double));
	// 	for (int j = 0; j < num; j++) {
	// 		A[i][j] = rand() % 10;
	// 		printf("%.2f ", A[i][j]);
	// 	}
	// 	printf("\n");
	// }

	// double *B = calloc(num, sizeof(double));

	// for (int i = 0; i < num; i++) {
	// 	B[i] = rand() % 100;
	// 	printf("%.2f ", B[i]);
	// }
	// printf("\n");


	// double t = wtime();
	// double *x = methodSimpleIteration(A, B, num); // простые итерации
	// t = wtime() - t;

	// printf("Time = %.5f\n", t);

	// printVector(x, num);

	double t = wtime();
	double *x = methodSeidels(A, B, num); // метод зейделя
	t = wtime() - t;

	printf("Time = %.5f\n", t);

	printVector(x, num);

	return 0;
}