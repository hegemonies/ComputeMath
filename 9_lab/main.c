#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define M_PI 3.14159265358979323846

// double Trigonometric_interpolation(double x_find, double *x, double *y, int n)
// {
// 	double complex res = 0.0 + 0.0 * I;

// 	double h = x[1] - x[0];

// 	// printf("%.3lf :: %.3lf\n", -((double)n/2) + 1, (double)n / 2);

// 	for (double j = -((double)n/2) + 1; j <= (double)n / 2; j++) {
// 	// for (double j = 0; j < n - 1; j++) {
// 		double complex tmp = 0.0 + 0.0 * I;
// 		for (double i = 0.0; i < n - 1; i++) {
// 			double complex y_ = y[(int)i];
// 			tmp += y_ * cexp((-2) * M_PI * I * ((i * j) / (double)n));
// 		}
// 		tmp /= n;

// 		res += tmp * cexp(2 * M_PI * I * j * ((x_find - x[0]) / ((double)n * h)));
// 		// res += tmp * cexp(2 * M_PI * I * j * x_find);
// 	}

// 	return creal(res);
// }

double Trigonometric_interpolation(double x_find, double *x, double *y, int n)
{
	double res = 0.0;
	double h = x[1] - x[0];

	for (int j = -(n/2) + 1, i = 0.0; j <= (n/2); j++, i++) {
		double tmp = 0.0;
		for (double k = 0.0; k < n; k++) {
			tmp += cos(M_PI * (-2) * (k * j / n)) * y[(int)k];
		}
		tmp /= n;

		res += cos(2 * M_PI * j *((x_find - x[0]) / (n * h))) * tmp;
	}

	return res;
}


int main()
{
	double y[3] = { 1.7321, 2.0, 2.2361 };
	double x[3] = { 3.0, 4.0, 5.0 };

	int n = 3;

	FILE *out = fopen("out.txt", "w");

	double k = 0;

	for (int i = 0; i < 40; i++) {
		fprintf(out, "%.3lf %.3lf\n", k, Trigonometric_interpolation(k, x, y, n));
		k += 0.2;
	}

	fclose(out);

	return 0;
}