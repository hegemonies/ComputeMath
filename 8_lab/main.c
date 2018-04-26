#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double Splines(double x_find, double *x, double *y, int n)
{
	double res = 0.0;

	int _i = 0;
	for (_i = 0; x_find > x[_i] && _i < n; _i++) { }

	double **c = calloc(0.0, sizeof(double*) * (n - 1));

	for (int i = 0; i < n - 1; i++) {
		c[i] = calloc(0.0, sizeof(double) * (n - 1));
		for (int j = 0; j < n - 1; j++) {
			if (i == j) {
				c[i][j] = ((x[i] - x[i - 1]) + (x[i + 1] - x[i])) / 3;
				continue;
			}
			if (j == i + 1) {
				c[i][j] = (x[i + 1] - x[i]) / 6;
				continue;
			}
			if (j == i - 1) {
				c[i][j] = (x[i] - x[i - 1]) / 6;
				continue;
			}
			c[i][j] = 0;
		}
	}

	double *d = calloc(0.0, sizeof(double) * (n - 1));

	for (int i = 1; i < n - 1; i++) {
		d[i] = ((y[i + 1] - y[i]) / (x[i + 1] - x[i])) - ((y[i] - y[i - 1]) / (x[i - 1] - x[i]));
	}

	double *M = calloc(0.0, sizeof(double) * (n - 1));

	for (int i = 0; i < n - 2; i++) {
		double tmp = 0.0;
		for (int j = ((i == 0) ? i : i - 1); j < ((i == n - 2) ? (i + 2) : (i + 3)); j++) {
			tmp += c[i][j];
		}
		M[i + 1] = d[i] / tmp;
	}


	res = M[_i - 1] * (pow(x[_i] - x_find, 3) / ((x[_i] - x[_i - 1]) * 6));
	res += M[_i]  * (pow(x_find - x[_i - 1], 3) / (6 * (x[_i] - x[_i - 1])));
	res += (y[_i - 1] - (M[_i - 1] * pow(x[_i] - x[_i - 1], 2)) / 6) * ((x[_i] - x_find) / (x[_i] - x[_i - 1]));
	res += (y[_i] - ((M[_i] * pow(x[_i] - x[_i - 1], 2)) / 6)) * ((x_find - x[_i - 1]) / (x[_i] - x[_i - 1]));

	for (int i = 0; i < n - 1; i++)
		free(c[i]);
	free(c);
	free(d);
	free(M);

	return res;
}

int main()
{
	double y[3] = { 1.7321, 2.0, 2.2361 };
	double x[3] = { 3.0, 4.0, 5.0 };

	int n = 3;
	// printf("%.4lf\n", Splines(3.6, x, y, n));

	FILE *out = fopen("out.txt", "w");

	double k = 0;

	for (int i = 0; i < 40; i++) {
		fprintf(out, "%.3lf %.3lf\n", k, Splines(k, x, y, n));
		k += 0.2;
	}

	fclose(out);

	return 0;
}