#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double eps = 1e-5;

double outD1 = 0.0;

double diff(double x, double y, double D1, double D2)
{
	if (x == 0) {
		x = 0.00001;
	}
	return pow(D2, 5) - cos(x) * D2 - sin(x) - 5 * log(x) * D1 - y * (x + 3);
}

double f(double x, double y, double D1)
{
	double a = 0, b = 2;
	double fa = 0, fb = 0;

	do {
		fa = diff(x, y, D1, a--);
		fb = diff(x, y, D1, b++);
	} while (fa * fb > 0);
	
	double c = 0;

	while(fabs(b - a) > eps) {
		c = (a + b) / 2;
		if(diff(x, y, D1, a) * diff(x, y, D1, c) < 0)
			b = c;
		else if (diff(x, y, D1, c) * diff(x, y, D1, b) < 0)
			a = c;
	}
	return (a + b) / 2;
}

double Runge_Kutt(double a, double b, double h, double y, double D1)
{
	double x = a;
	double y_t = 0.0;
	double D1_t = 0.0;

	while (x < b) {
		y_t = y + (h/2) * D1;
		D1_t = D1 + (h/2) * f(x, y, D1);
		y += h * D1_t;
		D1 += h * f(x + h/2, y_t, D1_t);
		x += h;
	}

	h = b - x;
	y_t = y + (h/2) * D1;
	D1_t = D1 + (h/2) * f(b, y, D1);
	y += h * D1_t;
	outD1 = D1 + h * f(b + h/2, y_t, D1_t);

	return y;
}

double MethodShooting(double x0, double x1, double y0, double y1, double h)
{
	double al = 1.0;
	double bt = 0.0;
	double fa = 0.0;
	double fb = 0.0;

	do {
		fa = Runge_Kutt(x0, x1, h, y0, al) - y1;
		fb = Runge_Kutt(x0, x1, h, y0, bt) - y1;
		al -= h;
		bt += h;
	} while (fa * fb > 0);

	double c = 0.0;
	while (fabs(bt - al) > eps) {
		c = (al + bt) / 2;
		if ((((Runge_Kutt(x0, x1, h, y0, al) - y1) * (Runge_Kutt(x0, x1, h, y0, c) - y1)) < 0)) {
			bt = c;
		} else if ((((Runge_Kutt(x0, x1, h, y0, c) - y1) * (Runge_Kutt(x0, x1, h, y0, bt) - y1)) < 0)) {
			al = c;
		}
	}

	return (al + bt) / 2;
}

// double Splines(double x_find, double *x, double *y, int n)
// {
// 	double res = 0.0;

// 	int _i = 0;
// 	for (_i = 0; x_find > x[_i] && _i < n; _i++) { }

// 	double **c = calloc(0.0, sizeof(double*) * (n - 1));

// 	for (int i = 0; i < n - 1; i++) {
// 		c[i] = calloc(0.0, sizeof(double) * (n - 1));
// 		for (int j = 0; j < n - 1; j++) {
// 			if (i == j) {
// 				c[i][j] = ((x[i] - x[i - 1]) + (x[i + 1] - x[i])) / 3;
// 				continue;
// 			}
// 			if (j == i + 1) {
// 				c[i][j] = (x[i + 1] - x[i]) / 6;
// 				continue;
// 			}
// 			if (j == i - 1) {
// 				c[i][j] = (x[i] - x[i - 1]) / 6;
// 				continue;
// 			}
// 			c[i][j] = 0;
// 		}
// 	}

// 	double *d = calloc(0.0, sizeof(double) * (n - 1));

// 	for (int i = 1; i < n - 1; i++) {
// 		d[i] = ((y[i + 1] - y[i]) / (x[i + 1] - x[i])) - ((y[i] - y[i - 1]) / (x[i - 1] - x[i]));
// 	}

// 	double *M = calloc(0.0, sizeof(double) * (n - 1));

// 	for (int i = 0; i < n - 2; i++) {
// 		double tmp = 0.0;
// 		for (int j = ((i == 0) ? i : i - 1); j < ((i == n - 2) ? (i + 2) : (i + 3)); j++) {
// 			tmp += c[i][j];
// 		}
// 		M[i + 1] = d[i] / tmp;
// 	}


// 	res = M[_i - 1] * (pow(x[_i] - x_find, 3) / ((x[_i] - x[_i - 1]) * 6));
// 	res += M[_i]  * (pow(x_find - x[_i - 1], 3) / (6 * (x[_i] - x[_i - 1])));
// 	res += (y[_i - 1] - (M[_i - 1] * pow(x[_i] - x[_i - 1], 2)) / 6) * ((x[_i] - x_find) / (x[_i] - x[_i - 1]));
// 	res += (y[_i] - ((M[_i] * pow(x[_i] - x[_i - 1], 2)) / 6)) * ((x_find - x[_i - 1]) / (x[_i] - x[_i - 1]));

// 	for (int i = 0; i < n - 1; i++)
// 		free(c[i]);
// 	free(c);
// 	free(d);
// 	free(M);

// 	return res;
// }


void set_h(double *h, double *X, int n)
{
	for (int i = 1; i < n; i++)
		h[i] = X[i] - X[i - 1];
}

void set_d(double *d, double *h, double *Y, int n)
{
	for (int i = 1; i < n - 1; i++)
		d[i] = ((Y[i + 1] - Y[i]) / h[i + 1]) - ((Y[i] - Y[i - 1]) / h[i]);
}

void set_C(double *C, double *h, int n)
{
	for (int i = 1; i < n - 1; i++) {
		for (int j = 1; j < n - 1; j++) {
			if (i == j) {
				C[i * n + j] = (h[i] + h[i + 1]) / 3;
			} else if (j == i + 1) {
				C[i * n + j] = h[i + 1] / 6;
			} else if (j == i - 1) {
				C[i * n + j] = h[i] / 6;
			} else {
				C[i * n + j] = 0;
			}
		}
	}
}

int Matrix_max_first_elem(double *a, int j, int n)
{
	double num = 0;
	int num_ind = 0;

	for (int z = j; z < n; z++) {
		if (fabsf(a[z * (n + 1) + j]) > num) {
			num = a[z * (n + 1) + j];
			num_ind = z;
		}
	}

	return num_ind;
}

void Matrix_swap_line(double *a, int j, int k, int n)
{
	double buf;
	for (int i = 1; i < n + 1; i++) {
		buf = a[j * (n + 1) + i];
		a[j * (n + 1) + i] = a[k * (n + 1) + i];
		a[k * (n + 1) + i] = buf;
	}
}

void Matrix_answer(double *arr_arg, double *a, int n)
{
	for (int i = n - 1; i > 0; i--) {
		for (int j = n - 1; j != i; j--) {
			a[i * (n + 1) + n] = a[i * (n + 1) + n] - (a[i * (n + 1) + j] * arr_arg[j]);
			a[i * (n + 1) + j] = 0;
		}

		arr_arg[i] = a[i * (n + 1) + n] / a[i * (n + 1) + i];
	}
}

void set_M(double *M, double *C, double *d, int n)
{
	double *arr = malloc(sizeof(double) * n * (n + 1));

	for (int i = 1; i < n - 1; i++) {
		for (int j = 1; j < n; j++) {
			arr[i * n + j] = C[i * n + j];
		}
		arr[i * n + (n - 1)] = d[i];
	}

	double mult;

	for (int j = 1; j < n - 2; j++) {

		int num_ind = Matrix_max_first_elem(arr, j, n - 1);
		Matrix_swap_line(arr, j, num_ind, n - 1);

		for (int i = j + 1; i < (n - 1); i++) {
			if (arr[i * (n - 1) + j] != 0) {
				mult = -(arr[i * (n - 1) + j] / arr[j * (n - 1) + i]);
			} else {
				break;
			}

			for (int k = j; k < n; k++) {
				arr[i * (n - 1) + k] += mult * arr[j * (n - 1) + k];
			}
		}
	}

	Matrix_answer(M, arr, (n - 1));
}

void print_M(double *M, int n)
{
	for (int i = 0; i < n; i++)
		printf("M[%d] = %.3f\n", i, M[i]);
}

int set_i(double *X, double x, int n)
{
	int i = 0;
	for (int j = 1; j < n; j++) {
		if (X[j - 1] <= x && x <= X[j])
			i = j;
	}

	return i;
}

double set_s(double *X, double *Y, double x, double *h, double *M, int i)
{
	double s = M[i - 1] * (pow((X[i] - x), 3) / (6 * h[i]));
	s += M[i] * (pow((x - X[i - 1]), 3) / (6 * h[i]));
	s += (Y[i - 1] - ((M[i - 1] * pow(h[i], 2)) / 6)) * ((X[i] - x) / h[i]);
	s += (Y[i] - ((M[i] * pow(h[i], 2)) / 6)) * ((x - X[i - 1]) / h[i]);

	return s;
}

double Splines(double *X, double *Y, double x, int n)
{
	double *h = malloc(sizeof(double) * n);
	double *d = malloc(sizeof(double) * (n - 1));
	double *C = malloc(sizeof(double) * n * n);
	double *M = malloc(sizeof(double) * n);

	set_h(h, X, n);

	set_d(d, h, Y, n);

	set_C(C, h, n);

	set_M(M, C, d, n);

	int i = set_i(X, x, n);

	double s = set_s(X, Y, x, h, M, i);

	free(M);
	free(C);
	free(d);
	free(h);

	return s;
}

double F(double x, double *y)
{
	if (x == 0)
		return y[0];
	if (x == 0.2)
		return y[1];
	if (x == 0.4)
		return y[2];
	if (x == 0.6)
		return y[3];
	if (x == 0.8)
		return y[4];
	if (x == 1.0)
		return y[5];
	return 0;
}

double Form_of_Simpson(double a, double b, double h, double *y)
{
	double res = 0.0;

	double j = a;

	for (int i = 1; j <= b - h; i++, j += h) {
		res += (i % 2 ? 4 : 2) * F(j, y);
	}

	res += F(a, y) + F(b, y);
	res = (res * h) / 3;

	return res;
}

double double_counting(double (*method)(double, double, double, double *), double a, double b, double h, double Eps, double *y)
{
	h = b - a;
	double prev = method(a, b, h, y);
	h /= 2;
	double cur = method(a, b, h, y);

	int count = 0;
	while (fabs(prev - cur) > Eps) {
		prev = cur;
		h /= 2;
		cur = method(a, b, h, y);

		count++;
	}
	printf("Count iteration = %d\n", count);

	return cur;
}

int main()
{
	double a = 0.0;
	double b = 1.0;
	double h = 0.2;

	int size = (a + b) / h;

	double y[size + 1];

	double x0 = 0.0;
	double y0 = 3.0;
	double x1 = 1.0;
	double y1 = 2.0;

	double D1 = MethodShooting(x0, x1, y0, y1, h);

	FILE *out = fopen("Runge_Kutt.txt", "w");
	
	int j = 0;
	for (double i = a; i <= b; i += h) {
		fprintf(out, "%.1lf %.5lf\n", i, Runge_Kutt(x0, i, h, y0, D1));
		y[j] = Runge_Kutt(x0, i, h, y0, D1);
		printf("%.1lf %.3lf\n", i, y[j]);
		j++;
	}
	printf("\n");

	fclose(out);

	double x[6] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };

	FILE *splines_out = fopen("Splines.txt", "w");
	
	for (double i = a; i <= b; i += 0.1) {
		double tmp = Splines(x, y, i, 6);
		fprintf(splines_out, "%.1lf %.5lf\n", i, tmp);
		printf("%.1lf %.3lf\n", i, tmp);
	}
	printf("\n");

	fclose(splines_out);

	double Eps = 1e-5;
	printf("Integral = %.8lf\n", double_counting(&Form_of_Simpson, a, b, h, Eps, y));
	// printf("My = %.10lf\n", Form_of_Trapeziums(a, b, h, y));
	// printf("Sanka =    %.10lf\n", SimpsonIntegr(a, b, y));

	return 0;
}