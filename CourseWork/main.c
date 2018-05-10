#include <stdio.h>
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

int main()
{
	double a = 0.0;
	double b = 1.0;
	double h = 0.2;

	double x0 = 0.0;
	double y0 = 3.0;
	double x1 = 1.0;
	double y1 = 2.0;

	int size = (a + b) / h;

	double y[size];
	y[0] = y0;
	y[size] = y1;

	double D1 = MethodShooting(x0, x1, y0, y1, h);

	FILE *out = fopen("Runge_Kutt.txt", "w");
	
	int j = 0;
	for (double i = a; i <= b; i += h) {
		fprintf(out, "%.1lf %.5lf\n", i, Runge_Kutt(x0, i, h, y0, D1));
		y[j] = Runge_Kutt(x0, i, h, y0, D1);
		printf("%.3lf\n", y[j]);
		j++;
	}

	fclose(out);

	return 0;
}