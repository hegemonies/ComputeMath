#include <stdio.h>
#include <stdlib.h>

// double P_Newton(double x_find, double *x, double *y, int n)
// {
// 	double res = 0.0;
// 	double res_mult;
// 	double res_add;
// 	for (int i = 1; i < n - 1; i++) {
// 		res_add = 0.0;
// 		for (int j = 0; j < i + 1; j++) {
// 			res_mult = 1.0;
// 			for (int k = 0; k < i + 1; k++) {
// 				if (j != k) {
// 					res_mult *= x[j] - x[k];
// 				}
// 			}
// 			res_add += y[j] / res_mult;
// 		}
// 		res_mult = 1.0;
// 		for (int j = 0; j < i; j++) {
// 			res_add *= (x_find - x[j]);
// 		}

// 		res += res_add;
// 	}


// 	return res + y[0];
// }

void make_diff_matx(double *res, double *y, int n)
{
	int size_line = n - 1;
	for (int i = 0; i < n - 1; i++) {
		for (int j = 0; j < (n - i - 1); j++) {
			if (i == 0) {
				res[j] = y[j + 1] - y[j];
			} else {
				res[i * size_line + j] = res[(i - 1) * size_line + (j + 1)] - res[(i - 1) * size_line + j];
			}
		}
	}
}

int factorial_iter(int p, int c, int max)
{
    if (c > max)
        return p;
    return factorial_iter(p * c, c + 1, max);
}

int factorial(int n)
{
    return factorial_iter(1, 1, n);
}

double P_Newton(double x_find, double *x, double *y, int n)
{
	double res = 0.0;
	double h = x[1] - x[0];

	double *diff_matx = calloc(0, sizeof(double) * (n - 1) * (n - 1));

	make_diff_matx(diff_matx, y, n);

	double q = (x_find - x[0]) / h;
	double tmp_q;

	for (int i = 1; i < n; i++) {
		tmp_q = 1.0;
		for (int j = 0; j < i; j++) {
			tmp_q *= q - j;
		}

		res += (diff_matx[(i - 1) * (n - 1)] / factorial(i)) * tmp_q;
	}

	return res + y[0];
}

int main()
{
	// double y[6] = { -1, 2, 13, 44, 107, 214 };
	// int n = 6;
	// double *res = malloc(sizeof(double) * 5 * 5);

	// for (int i = 0; i < n - 1; i++) {
	// 	for (int j = 0; j < (n - i - 1); j++) {
	// 		res[i * (n - 1) + j] = 0.0;
	// 	}
	// }

	// make_diff_matx(res, y, n);

	// for (int i = 0; i < n - 1; i++) {
	// 	for (int j = 0; j < (n - i - 1); j++) {
	// 		printf("%.3lf\n", res[i * (n - 1) + j]);
	// 	}
	// 	printf("\n");
	// }

	// for (int i = 0; i < n - 1; i++)
	// 	free(res[i]);
	// free(res);
	
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