#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f(double x)
{
	return sqrt(x);
}

double df(double x)
{
	return (1 / (2 * sqrt(x)));
}

double df2(double x)
{
	return (3 / (8 * (pow(x, 5 / 2))));
}

double df3(double x)
{
	return (105 / (32 * (pow(x, 9 / 2))));
}

double P(double x_find, double *x, double *y, int n)
{
	if (!x) {
		return -1;
	}

	double *res = calloc(n, sizeof(double));
	double res_product = 0.0;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			res_product *= (x_find - x[j]) / (x[i] - x[j]);
		}
		res[i] = res_product * y[i];
	}

	res_product = 0.0;

	for (int i = 0; i < n; i++) {
		res_product += res[i];
	}

	free(res);

	return res_product;
}

// double outputInFile(double x, double y) // TODO
// {

// }

int main()
{
	double y[3] = { 1.7321, 2.0, 2.2361 };
	double x[3] = { 3.0, 4.0, 5.0 };

	int n = 3;

	printf("%lf\n", P(4.41, x, y, n));

	return 0;
}