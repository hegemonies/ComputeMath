#include <stdio.h>
#include <stdlib.h>

double P_Newton(double x_find, double *x, double *y, int n)
{
	double res = 0.0;
	double res_mult;
	double res_add;
	for (int i = 1; i < n - 1; i++) {
		res_add = 0.0;
		for (int j = 0; j < i + 1; j++) {
			res_mult = 1.0;
			for (int k = 0; k < i + 1; k++) {
				if (j != k) {
					res_mult *= x[j] - x[k];
				}
			}
			res_add += y[j] / res_mult;
		}
		res_mult = 1.0;
		for (int j = 0; j < i; j++) {
			res_add *= (x_find - x[j]);
		}

		res += res_add;
	}


	return res + y[0];
}

int main()
{
	double y[3] = { 1.7321, 2.0, 2.2361 };
	double x[3] = { 3.0, 4.0, 5.0 };

	int n = 3;

	FILE *out = fopen("out.txt", "w");

	double k = 0;

	for (int i = 0; i < 40; i++) {
		fprintf(out, "%.3lf %.3lf\n", k, P_Newton(k, x, y, n));
		k += 0.2;
	}

	fclose(out);

	return 0;
}