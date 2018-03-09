#include <stdio.h>
#include <math.h>

#define eps 0.001

double f(double x) // y = 4x^2 - 4x + 1
{
	return x*x - 3;
}

double df(double x) // y' = 8x + 4
{
	return 2*x;
}

double hdm(double a, double b)
{
	double c = (a + b) / 2;

	while (fabs((b - a) / 2) > eps) {
		if (f(a) * f(c) < 0) {
			b = c;
		} else {
			a = c;
		}

		c = (a + b) / 2;
		printf("che1\n");
	}

	return c;
}

double mc(double a, double b)
{
	double tmp = 99999;
	double c = (a*f(b) - b*f(a)) / (f(b) - f(a));

	while (fabs(tmp - c) > eps) {
		if (f(a) * f(c) < 0) {
			b = c;
		} else {
			a = c;
		}

		tmp = c;
		c = (a*f(b) - b*f(a)) / (f(b) - f(a));
		printf("che2\n");
	}

	return c;
}

double mn(double a, double b)
{
	double h = 0;
	double x = 1;

	do {
		h = f(x) / df(x);
		x = x - h;
		printf("che3\n");
	} while (fabs(h) > eps);

	return x;
}

int main()
{
	double a = -1; // начало заданного отрезка
	double b = 3; // конец

	printf("f(a) * f(b) = %.3lf\n\n", df(a) * df(b));

	double c = hdm(a, b);

	printf("x = %.4lf\n", c);
	printf("f(x) = %.4lf\n\n", f(c));

	c = mc(a, b);

	printf("x = %.4lf\n", c);
	printf("f(x) = %.4lf\n\n", f(c));

	c = mn(a, b);

	printf("x = %.4lf\n", c);
	printf("f(x) = %.4lf\n", f(c));

	return 0;
}