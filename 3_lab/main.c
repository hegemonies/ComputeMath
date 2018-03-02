#include <stdio.h>
#include <math.h>

#define eps 0.1

double f(double x) // y = 4x^2 - 4x + 1
{
	return (4*x*x) - (4*x) + 1;
}

double df(double x) // y' = 8x + 4
{
	return 8*x - 4;
}

double hdm(double a, double b)
{
	if (df(a) * df(b) > 0) {
		return -1;
	}

	double c = fabs(b - a) / 2;

	while (f(c) > eps) {

		if (f(a) * f(c) < 0) {
			b = c;
		} else {
			a = c;
		}
		c = fabs(b - a) / 2;
	}

	return c;
}

double mc(double a, double b)
{
	if (df(a) * df(b) > 0) {
		return -1;
	}

	double next = 0;
	double tmp = b;

	while (fabs(next - b) > eps) {
		tmp = next;
		next = a - (df(a) / (df(a) - df(tmp))) * (a - tmp);
		b = tmp;
	}

	return next;
}

double mn(double a, double b)
{
	double h = 0;
	double x = 1;

	do {
		h = f(x) / df(x);
		x = x - h;
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