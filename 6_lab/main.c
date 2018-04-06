#include <stdio.h>

double P(double x_find, double *x, double *y, int n0, int nn)
{
	if ((nn - n0) == 1) {
		return (y[n0] * (x_find - x[nn]) -y[nn] * (x_find - x[n0])) / (x[n0] - x[nn]);
	} else {
		double P1 = P(x_find, x, y, n0, nn - 1);
		double P2 = P(x_find, x, y, n0 + 1, nn);
		return (((P1 * (x_find - x[nn])) - (P2 * (x_find - x[n0]))) / (x[n0] - x[nn]));
	}
}

int main()
{
	double y[3] = { 1.7321, 2.0, 2.2361 };
	double x[3] = { 3.0, 4.0, 5.0 };

	int n = 3;

	FILE *out = fopen("out.txt", "w");

	double k = 0.0;

	for (int i = 0; i < 40; i++) {
		fprintf(out, "%.3lf %.3lf\n", k, P(k, x, y, 0, n));
		k += 0.2;
	}

	fclose(out);

	return 0;
}