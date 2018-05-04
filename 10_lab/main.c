#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 1000000

double F(double x)
{
	return sqrt(x);
}

double Form_of_Trapeziums(double a, double b, double h)
{
	double res = 0.0;

	double j = a;
	for (int i = 1; j <= b; i++, j += h) {
		res += F(j);
	}

	res += (F(a) + F(b)) / 2;
	res *= h;

	return res;
}

double Form_of_Simpson(double a, double b, double h)
{
	double res = 0.0;

	double j = a;

	for (int i = 1; j <= b - h; i++, j += h) {
		res += (i % 2 ? 4 : 2) * F(j);
	}

	res += F(a) + F(b);
	res = (res * h) / 3;

	// res = ((b - a) / 6) * (F(a) + 4 * F((a + b) / 2) + F(b));

	return res;
}

double double_counting(double (*method)(double, double, double), double a, double b, double h, double eps)
{
	double prev = method(a, b, h);
	h /= 2;
	double cur = method(a, b, h);

	int count = 0;
	while (fabs(prev - cur) > eps) {
		prev = cur;
		h /= 2;
		cur = method(a, b, h);

		count++;
	}
	printf("Count iteration = %d\n", count);

	return cur;
}

int main()
{
	double a = 1.0;
	double b = 2.0;
	double h = (b - a) / N;

	printf("Result (Trapeziums) simple = %.8lf\n", Form_of_Trapeziums(a, b, h));
	printf("Result (Simpson) simple = %.8lf\n\n", Form_of_Simpson(a, b, h));

	double eps = 1e-4;
	printf("Result (double_counting Trap) for eps %lf = %.8lf\n", eps, double_counting(&Form_of_Trapeziums, a, b, h, eps));
	printf("Result (double_counting Simpson) for eps %lf = %.8lf\n\n", eps, double_counting(&Form_of_Simpson, a, b, h, eps));

	eps = 1e-5;
	printf("Result (double_counting Trap) for eps %.10lf = %.8lf\n", eps, double_counting(&Form_of_Trapeziums, a, b, h, eps));
	printf("Result (double_counting Simpson) for eps %.10lf = %.8lf\n\n", eps, double_counting(&Form_of_Simpson, a, b, h, eps));

	eps = 1e-6;
	printf("Result (double_counting Trap) for eps %.10lf = %.8lf\n", eps, double_counting(&Form_of_Trapeziums, a, b, h, eps));
	printf("Result (double_counting Simpson) for eps %.10lf = %.8lf\n\n", eps, double_counting(&Form_of_Simpson, a, b, h, eps));

	eps = 1e-7;
	printf("Result (double_counting Trap) for eps %.10lf = %.8lf\n", eps, double_counting(&Form_of_Trapeziums, a, b, h, eps));
	printf("Result (double_counting Simpson) for eps %.10lf = %.8lf\n\n", eps, double_counting(&Form_of_Simpson, a, b, h, eps));

	eps = 1e-8;
	printf("Result (double_counting Trap) for eps %.10lf = %.8lf\n", eps, double_counting(&Form_of_Trapeziums, a, b, h, eps));
	printf("Result (double_counting Simpson) for eps %.10lf = %.8lf\n\n", eps, double_counting(&Form_of_Simpson, a, b, h, eps));

	return 0;
}