#include <stdio.h>
#include <stdlib.h>

double Polinom_Newton(double x_find, double *x, double *y, int n)
{
	double res_add = 0.0;
	double res_mult = 1.0;
	// double a = 0.0;

	// for (int m = 0; m =< n; m++) { // Pn(x)
	// 	for (int i = 0; i < m; i++) { // Ai
	// 		for (int j = 0; j =< i; j++) {
	// 			res_mult *= x[i] - x[j];
	// 		}
	// 		res_add += y[i] / res_mult;
	// 	}
	// 	a = res_add * (x_find - x[m]);
	// }

	double res = 0.0;

	for (int m = 1; m < n - 1; m++) {
		for (int i = 1; i < m + 1; i++) {
			for (int j = 1; j < m + 1; j++) {
				if (i != j) {
					res_mult *= x[i] - x[j];
				}
			}
			res_add += y[i] / res_mult;
		}
		for (int i = 0; i < m; i++) {
			res = x_find - x[i];
		}
		res = res_add * res;
		// printf("check\n");
	}
	// printf("\n");
	

	// return res_add;
	return res + y[0];
}

int main()
{
	double y[3] = { 1.7321, 2.0, 2.2361 };
	double x[3] = { 3.0, 4.0, 5.0 };

	int n = 3;

	FILE *out = fopen("out.txt", "w");

	double k = 3.2;

	for (int i = 0; i < 10; i++) {
		fprintf(out, "%.3lf %.3lf\n", k, Polinom_Newton(k, x, y, n));
		k += 0.2;
	}

	fclose(out);

	return 0;
}