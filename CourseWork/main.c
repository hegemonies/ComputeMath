#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double eps = 1e-6;

double outD1 = 0.0;

double diff(double x, double y, double D1, double D2)
{
	if (x == 0) {
		x = 0.0000000000000000000000001;
	}
	return pow(D2, 5) - cos(x) * D2 - sin(x) - 5 * log(x) * D1 - y * (x + 3);
}

double *addition_of_vectors(double *v1, double *v2, int n)
{
	double *v3 = malloc(sizeof(double) * n);

	for (int i = 0; i < n; i++) {
		v3[i] = v1[i] + v2[i];
	}

	return v3;
}

double *multiple_dig_by_vector(double a, double *v, int n)
{
	double *v2 = malloc(sizeof(double) * n);

	for (int i = 0; i < n; i++) {
		v2[i] = a * v[i];
	}

	return v2;
}

double *f(double x, double *y)
{
	double a = 0, b = 2;
	double fa = 0, fb = 0;

	do {
		fa = diff(x, y[0], y[1], a--);
		fb = diff(x, y[0], y[1], b++);
	} while (fa * fb > 0);
	
	double c = 0;

	while(fabs(b - a) > eps) {
		c = (a + b) / 2;
		if(diff(x, y[0], y[1], a) * diff(x, y[0], y[1], c) < 0)
			b = c;
		else if (diff(x, y[0], y[1], c) * diff(x, y[0], y[1], b) < 0)
			a = c;
	}

	double *tmp_y = malloc(sizeof(double) * 2);
	tmp_y[0] = y[1];
	tmp_y[1] = (a + b) / 2;
	
	return  tmp_y;
}

double *Runge_Kutt(double a, double b, double h, double *y0)
{
	double *tmp_y;
	double *y = y0;
	for (double i = a + h; i <= b; i += h) {
		tmp_y = addition_of_vectors(y, multiple_dig_by_vector(h, f(i, y), 2), 2);
		y = addition_of_vectors(y, multiple_dig_by_vector(h / 2, addition_of_vectors(f(i, y), f(i + h, tmp_y), 2), 2), 2);
	}

	return y;
}


double MethodShooting(double x0, double x1, double y0, double y1, double h)
{
	double al = 1.0;
	double bt = 0.0;
	double fa = 0.0;
	double fb = 0.0;

	double tmp[2];
	double *vt;

	do {
		tmp[0] = y0;
		tmp[1] = al;
		vt = Runge_Kutt(x0, x1, h, tmp);
		fa = vt[0] - y1;
		tmp[1] = bt;
		vt = Runge_Kutt(x0, x1, h, tmp);
		fb = vt[0] - y1;
		al -= h;
		bt += h;
	} while (fa * fb > 0);

	double c = 0.0;
	while (fabs(bt - al) > eps) {
		c = (al + bt) / 2;

		tmp[1] = al;
		double *tmp1 = Runge_Kutt(x0, x1, h, tmp);
		tmp[1] = c;
		double *tmp2 = Runge_Kutt(x0, x1, h, tmp);
		tmp[1] = bt;
		double *tmp4 = Runge_Kutt(x0, x1, h, tmp);

		if ((((tmp1[0] - y1) * (tmp2[0] - y1)) < 0)) {
			bt = c;
		} else if ((((tmp2[0] - y1) * (tmp4[0] - y1)) < 0)) {
			al = c;
		}
	}

	return (al + bt) / 2;
}

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

double Form_of_Simpson(double a, double b, double h, double *y)
{
	double res = 0.0;
	int n = (int)((a + b) / h) + 1;

	for (int i = 1; i < 2*n - 1; i += 2)
		res += 4 * y[i];

	for (int i = 2; i < n - 2; i += 2) 
		res += 2 * y[i];

	res += y[0] + y[n];
	res = res * h / 3;


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


double **dbl_counting_Runge(double a, double b, double h, double Eps, double *y)
{
	int check = 1;
	double **y_prev;
	double **y_cur;
	do {
		int count_elem = (b - a) / h;
		printf("count_elem = %d\n", count_elem);
		y_prev = malloc(sizeof(double*) * count_elem);
		double t = a;
		for (int i = 0; t < b; i++, t += h) {
			y_prev[i] = malloc(sizeof(double*) * 2);
			y_prev[i] = Runge_Kutt(a, t, h, y);
		}

		h /= 2;
		count_elem = (b - a) / h;
		printf("count_elem = %d\n", count_elem);
		t = a;
		y_cur = malloc(sizeof(double*) * count_elem);
		for (int i = 0; t < b; i++, t += h) {
			y_cur[i] = malloc(sizeof(double*) * 2);
			y_cur[i] = Runge_Kutt(a, t, h, y);
		}

		double mod = 0.0;
		int j = 0;
		for (int i = 0; j < count_elem * 2; i++, j += 2) {
			mod = y_prev[i] - y_cur[j];
			if (fabs(mod) > Eps) {
				break;
			}
		}
		if (j == (count_elem * 2) - 1) {
			check = 0;
		}

		printf("chee\n");
	} while (check);

	return y_cur;
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
	printf("D1 = %.3lf\n", D1);

	FILE *out = fopen("Runge_Kutt.txt", "w");
	
	// double dy[size + 1];

	// int j = 0;
	printf("x\ty(x)\ty\'(x)\n");
	// for (double i = a; i <= b; i += h) {
	// 	double tmp[2] = { y0, D1 };
	// 	// double *tv = Runge_Kutt(x0, i, h, tmp);
	// 	double *tv = dbl_counting_Runge(x0, i, h, tmp);
	// 	fprintf(out, "%.1lf %.5lf\n", i, tv[0]);
	// 	y[j] = tv[0];
	// 	dy[j] = tv[1];
	// 	printf("%.1lf\t%.3lf\t%.3lf\n", i, y[j], dy[j]);
	// 	j++;
	// }
	double Eps = 1e-2;

	double tmp[2] = { y0, D1 };
	double **yt = dbl_counting_Runge(a, b, h, Eps, tmp);
	for (int i = 0; i < 6; i++) {
			printf("\t%.3lf\t%.3lf\n", yt[i][0], yt[i][1]);
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


	// printf("I = %.10lf\n", Form_of_Simpson(a, b, h, y));
	printf("I = %.10lf\n", double_counting(Form_of_Simpson, a, b, h, Eps, y));

	return 0;
}