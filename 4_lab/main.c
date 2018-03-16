#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>

#define Eps 0.00001

typedef double (*F)(double, double, double);
typedef double (*W)(double, double, double);

double f1(double x, double y, double z)
{
	return (x * x) + (y * y) + (z * z);
}

double f2(double x, double y, double z)
{
	return 2*(x * x) + (y * y) + 4*z;
}

double f3(double x, double y, double z)
{
	return 3*(x * x) - 4*y + (z * z);
}

double df1(double x, double y, double z)
{
	return 2*x + 2*y + 2*z;
}

double df2(double x, double y, double z)
{
	if (z != 0) {
		return -4;
	}
	return 4*x + 2*y;
}

double df3(double x, double y, double z)
{
	if (y != 0) {
		return -4;
	}
	return 6*x + 2*z;
}

double *compF(F *f, double *vect, double *b, int n)
{
	double *tmp = calloc(n, sizeof(double));

	for (int i = 0; i < n; i++) {
		tmp[i] = f[i](vect[0], vect[1], vect[2]) - b[i];
	}

	return tmp;
}

double **compW(W *f, double *vect, int n)
{
	double **tmp = calloc(n, sizeof(double*));

	int check = 1;

	for (int i = 0; i < n; i++) {
		tmp[i] = calloc(n, sizeof(double));
		for (int j = 0; j < n; j++) {
			tmp[i][j] = f[i]((check == 1) ? vect[0] : 0, (check == 2) ? vect[1] : 0, (check == 3) ? vect[2] : 0);
			check++;
		}
		check = 1;
	}

	return tmp;
}

void swap(double **arr, int orig, int dest, int n)
{
	double tmp;
	for (int i = 0; i < n; i++) {
		tmp = arr[orig][i];
		arr[orig][i] = arr[dest][i];
		arr[dest][i] = tmp;
	}
}


double minor(double **matx)
{
	return (matx[0][0] * matx[1][1]) - (matx[0][1] * matx[1][0]);
}

double determinant(double **a, double **A, int n)
{
	double sum = 0;

	for (int i = 0; i < n; i++) {
		sum += a[0][i] * A[i][0];
	}

	return sum;
}

void printMatx(double **m, int n)
{
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%.3f ", m[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

double **reversMatx(double **matx, int n)
{
	double **mRev = calloc(n, sizeof(double));
	for (int i = 0; i < n; i++) {
		mRev[i] = calloc(n, sizeof(double));
	}

	//double deter = determinant(matx, n);

	double **matxTmp2 = calloc(2, sizeof(double));
	matxTmp2[0] = calloc(2, sizeof(double));
	matxTmp2[1] = calloc(2, sizeof(double));

	int i2 = 0;
	int j2 = 0;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				for (int m = 0; m < n; m++) {
					if (k != i && m != j) {
						matxTmp2[i2][j2] = matx[k][m];
						j2++;
						if (j2 > 1) {
							j2 = 0;
							i2++;
						}
					}
				}
			}
			// printMatx(matxTmp2, 2);
			mRev[j][i] = pow(-1, (i + 1) + (j + 1)) * minor(matxTmp2);
			i2 = 0;
			j2 = 0;
		}
		//printf("\n");
	}
	// printf("\n");

	// printMatx(mRev, n);


	// printf("deter = %.2f\n", determinant(matx, mRev, n));
	double deter = determinant(matx, mRev, n);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			mRev[i][j] /= deter;
		}
	}

	return mRev;
}

double **multipMatxOnMatx(double **a, double **b, int n)
{
	double **c = calloc(n, sizeof(double));
	for (int i = 0; i < n; i++) {
		c = calloc(n, sizeof(double));
	}

	for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            for (int j = 0; j < n; j++) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }

	return c;
}

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

double *vectorMinusVector(double *l, double *r, int n)
{
	double *res = calloc(n, sizeof(double));

	for (int i = 0; i < n; i++) {
		res[i] -= l[i] - r[i];
	}

	return res;
}

int main()
{
	int n = 3;

	F f[3];

	f[0] = f1;
	f[1] = f2;
	f[2] = f3;

	double B[3] = { 1, 0, 0 };

	double X0[3] = { 0.5, 0.5, 0.5 };

	W w[3];

	w[0] = df1;
	w[1] = df2;
	w[2] = df3;

	double *tmp = calloc(n, sizeof(double));
	double *prev = calloc(n, sizeof(double));
	double *x = calloc(n, sizeof(double));

	double *y = calloc(n, sizeof(double));

	double **tmp_yakobi;

	double *fxk;

	// int count = 0;

	double max = INT_MIN;
	x = X0;

	do {
		tmp_yakobi = compW(w, x, n);
		fxk = compF(f, x, B, n);
		y = matrxMultVect(reversMatx(tmp_yakobi, n), fxk, n);

		tmp = vectorMinusVector(x, prev, n);

		for (int i = 0; i < n; i++) {
			if (tmp[i] >= max) {
				printf("tmp[%d] = %.3f\n", i, tmp[i]);
				printf("max = %.3f\n", max);
				max = fabs(tmp[i]);
			}
		}

		x = vectorMinusVector(x, y, n);
		// printf("che = %d\n", max);
	} while (fabs(max) > Eps);

	printf("%.3f\n", x[0]);
	printf("%.3f\n", x[1]);
	printf("%.3f\n", x[2]);


	// double **tst = calloc(n, sizeof(double));
	// for (int i = 0; i < n; i++) {
	// 	tst[i] = calloc(n, sizeof(double));
	// 	for (int j = 0; j < n; j++) {
	// 	}
	// }

	// tst[0][0] = 2;
	// tst[0][1] = 1;
	// tst[0][2] = -1;

	// tst[1][0] = 3;
	// tst[1][1] = 2;
	// tst[1][2] = -2;

	// tst[2][0] = 1;
	// tst[2][1] = -1;
	// tst[2][2] = 2;

	// printf("det = %.2f\n", determinant(tst, n));

	// for (int i = 0; i < n; i++) {
	// 	for (int j = 0; j < n; j++) {
	// 		printf("%.2f ", tst[i][j]);
	// 	}
	// 	printf("\n");
	// }
	// printf("\n");

	// reversMatx(tst, n);

	// double **m = calloc(2, sizeof(double));
	// m[0] = calloc(2, sizeof(double));
	// m[1] = calloc(2, sizeof(double));

	// m[0][0] = 1;
	// m[0][1] = 1;
	// m[1][0] = 4;
	// m[1][1] = 5;

	// printf("deter = %.2f\n", deterFor2(m));

	return 0;
}