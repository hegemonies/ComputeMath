#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double P(double x_find, double *x, double *y, int n)
{
	if (!x) {
		return -1;
	}

	double *res = calloc(n, sizeof(double));
	double res_product;

	for (int i = 0; i < n; i++) {
		res_product = 1.0;
		for (int j = 0; j < n; j++) {
			if (i != j) {
				res_product *= (x_find - x[j]) / (x[i] - x[j]);
			}
		}
		res[i] = res_product * y[i];
		// printf("%f\n", res[i]);
	}

	res_product = 0.0;

	for (int i = 0; i < n; i++) {
		res_product += res[i];
	}

	free(res);

	return res_product;
}

int main()
{
	double y[3] = { 1.7321, 2.0, 2.2361 };
	double x[3] = { 3.0, 4.0, 5.0 };

	int n = 3;

	FILE *out = fopen("out.txt", "w");

	double k = 0.0;

	for (int i = 0; i < 40; i++) {
		fprintf(out, "%.3lf %.3lf\n", k, P(k, x, y, n));
		k += 0.2;
	}

	fclose(out);

	return 0;
}