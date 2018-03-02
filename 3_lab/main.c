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
	double c = fabs(b - a) / 2;

	while (f(c) > eps) {
		printf("a = %.4lf\n", a);
		printf("b = %.4lf\n", b);
		printf("c = %.4lf\n\n", c);

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
	double next = 0;
	double tmp = b;

	while (fabs(next - b) > eps) {
		tmp = next;
		// next = b - f(b) * (a - b) / (f(a) - f(b));
		next = a - (f(a) / (f(a) - f(tmp))) * (a - tmp);
		// a = b;
		b = tmp;
		printf("next = %.3lf\n", next);
		printf("a = %.3lf\n", a);
		printf("b = %.3lf\n\n", b);
	}

	return next;
}

int main()
{
	double a = -1; // начало заданного отрезка
	double b = 3; // конец

	printf("f(a) * f(b) = %.3lf\n", df(a) * df(b));

	double c = hdm(a, b);

	printf("x = %.4lf\n", c);
	printf("f(x) = %.4lf\n", f(c));

	c = mc(a, b);

	printf("x = %.4lf\n", c);
	printf("f(x) = %.4lf\n", f(c));

	return 0;
}