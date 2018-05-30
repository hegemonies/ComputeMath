#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>

using namespace std;

double eps = 1e-4;

double diff(double x, double y, double D1, double D2)
{
	if (x == 0) {
		x = 0.0001;
	}
	return pow(D2, 5) - cos(x) * D2 - sin(x) - 5 * log(x) * D1 - y * (x + 3);
}

double *addition_of_vectors(double *v1, double *v2, int n)
{
	double *v3 = new double[n];

	for (int i = 0; i < n; i++) {
		v3[i] = v1[i] + v2[i];
	}

	return v3;
}

double *multiple_dig_by_vector(double a, double *v, int n)
{
	double *v2 = new double[n];

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

	while(fabs(b - a) >= eps) {
		c = (a + b) / 2;
		if(diff(x, y[0], y[1], a) * diff(x, y[0], y[1], c) < 0)
			b = c;
		else if (diff(x, y[0], y[1], c) * diff(x, y[0], y[1], b) < 0)
			a = c;
	}

	double *tmp_y = new double[2];
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

double **dbl_counting_Runge(double a, double b, double h, double Eps, double *y, int *count, double *h_)
{
	double **y_prev;
	double **y_cur;
	double max;
	int count_elem;
	do {
		count_elem = (b - a) / h;
		y_prev = new double*[count_elem];
		double t = a;
		for (int i = 0; fabs(t - b) >= 1e-8; i++, t += h) {
			y_prev[i] = new double[2];
			y_prev[i] = Runge_Kutt(a, t, h, y);
		}

		h /= 2;
		count_elem = (b - a) / h;
		t = a;
		y_cur = new double*[count_elem];
		for (int i = 0; fabs(t - b) >= 1e-8; i++, t += h) {
			y_cur[i] = new double[2];
			y_cur[i] = Runge_Kutt(a, t, h, y);
		}

		double mod = 0.0;
		max = 0;
		int j = 0;
		int del = count_elem;
		int i;
		for (i = 0; j < del + (del % 2 ? 1 : 0); i++, j += 2) {
			mod = fabs(y_prev[i][0] - y_cur[j][0]);
			if (mod > max) {
				max = mod;
			}
		}
	} while (fabs(max) >= Eps);

	if (count) {
		*count = count_elem;
	}

	if (h_) {
		*h_ = h;
	}

	for (int i = 0, j = 0; i < count_elem; i += 2, j++) {
		y_cur[i][0] -= (y_cur[i][0] - y_prev[j][0]) / 3;
		y_cur[i][1] -= (y_cur[i][1] - y_prev[j][1]) / 3;
	}

	return y_cur;
}

double MethodShooting(double x0, double x1, double y0, double y1, double h)
{
	double al = 1.0;
	double bt = 0.0;
	double fa = 0.0;
	double fb = 0.0;

	double tmp[2];
	tmp[0] = y0;
	double *vt;
	do {
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
	double *tmp1;
	double *tmp2;
	double *tmp4;
	do {

		tmp[1] = al;
		tmp1 = Runge_Kutt(x0, x1, h, tmp);
		tmp[1] = c;
		tmp2 = Runge_Kutt(x0, x1, h, tmp);
		tmp[1] = bt;
		tmp4 = Runge_Kutt(x0, x1, h, tmp);

		if ((((tmp1[0] - y1) * (tmp2[0] - y1)) < 0)) {
			bt = c;
		} else if ((((tmp2[0] - y1) * (tmp4[0] - y1)) < 0)) {
			al = c;
		}

		c = (al + bt) / 2;
	} while (fabs(y1 - (tmp1[0] + tmp2[0] + tmp4[0]) / 3) > 1e-4);

	return ((al + bt) / 2);
}

double Form_of_Simpson(double a, double b, double h, double *y0)
{
	double res = 0.0;
	double j = a;

	double *tmp;

	for (int i = 1; j <= b - h; i++, j += h) {
		tmp = Runge_Kutt(a, j, h, y0);
		res += (i % 2 ? 4 : 2) * tmp[0];
	}

	tmp = Runge_Kutt(a, 0, h, y0);
	res += tmp[0];
	tmp = Runge_Kutt(a, b, h, y0);
	res += tmp[0];
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
	// printf("h = %lf\n", h);
	// printf("Count iteration = %d\n", count);

	return cur;
}


double dbl_count_D1(double x0, double x1, double y0, double y1, double h, double Eps)
{
	double prev;
	double cur;
	do {
		prev = MethodShooting(x0, x1, y0, y1, h);
		h /= 2;
		cur = MethodShooting(x0, x1, y0, y1, h);
	} while (fabs(prev - cur) >= Eps);
	
	return cur;
}

double h_i(int i, const pair<double, double> * xy)
{
	return xy[i].first - xy[i-1].first;
}

double b_i(int i, const pair<double, double> * xy)
{
	return (h_i(i, xy) + h_i(i+1, xy)) / 3;
}

double g_i(int i, const pair<double, double> * xy)
{
	return h_i(i, xy) / 6;
}

double d_i(int i, const pair<double, double> * xy)
{
	return ((xy[i+1].second - xy[i].second) / h_i(i+1, xy) - 
				(xy[i].second - xy[i-1].second) / h_i(i, xy));
}

vector <double> thomas_come_on(vector <double>  a, vector <double> b, 
									vector <double>  g, vector <double>  d, int n)
{
	vector <double> c(n);
	vector <double> p(n);
	vector <double> q(n);

	for(int i = 0; i < n; i++) {
		if (i == 0) {
			p[i] = -g[i] / b[i];
			q[i] =  d[i] / b[i];
			continue;
		}
		p[i] = g[i] / (-b[i] - a[i] * p[i-1]);
		q[i] = (a[i] * q[i-1] - d[i])/(-b[i] - a[i] * p[i-1]);
	}

	for(int i = n - 1 ; i >= 0; i--) {
		if( i == n) {
			c[i] = (a[i] * q[i-1] - d[i]) / (-b[i] - a[i] * p[i-1]);
			continue;
		}
		c[i] = p[i] * c[i + 1] + q[i];
	}
	return c;
}

int binary_search (pair<double, double> * xy, double x, int n)
{
	int idx = 0;
	if(x <= xy[0].first) {
		idx = 1;
	}
	else if (x >= xy[n-1].first) {
		idx = n-1;
	}
	else {
		int i = 0, j = n-1;
		while(i + 1 < j) {
			int k = i + (j - i) / 2;
			if (x <= xy[k].first) {
				j = k;
			} else {
				i = k;
			}
		}
		idx = j;
	}
	return idx;
}

double spline_eval(int i, double * M, pair<double, double> * v, double x)
{
	double s1 = M[i-1] * pow(v[i].first - x, 3) / (6 * h_i(i, v));

	double s2 = M[i] * pow(x - v[i-1].first, 3) / (6 * h_i(i, v));

	double s3 = (v[i-1].second - M[i-1] * pow(h_i(i, v), 2) / 6) * 
			(v[i].first - x) / h_i(i, v);

	double s4 = (v[i].second - M[i] * pow(h_i(i, v), 2) / 6) *
			(x - v[i-1].first) / h_i(i, v);

	return s1 + s2 + s3 + s4;
}

double cubic(double x, int n, pair<double, double> * v)
{
	int i = 0;
	double s = 0;
	vector <double> M(n - 2);
	vector <double> a(n - 2), b(n - 2), g(n - 2), d(n - 2);

	a[0] = g[n-3] = 0;
	for(int i = 1; i < n - 1; i++)
		b[i - 1] = b_i(i, v);
	for(int i = 2; i < n - 1; i++)
		a[i - 1] = g_i(i, v);
	for(int i = 0; i < n - 3; i++)
		g[i] = g_i(i+1, v);
	for(int i = 1; i < n - 1; i++)
		d[i-1] = d_i(i, v);

	M = thomas_come_on(a, b, g, d, n-2);
	i = binary_search(v, x, n);
	s = spline_eval(i, M.data(), v, x);

	return s;
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

	double D1 = dbl_count_D1(x0, x1, y0, y1, h, 1e-3);
	printf("D1 = %.3lf\n", D1);

	FILE *out = fopen("Runge_Kutt.txt", "w");
	printf("x\ty(x)\ty\'(x)\n");
	double Eps = 1e-3;

	double tmp[2] = { y0, D1 };
	int count_elem;
	double h_;
	double **yt = dbl_counting_Runge(a, b, h, Eps, tmp, &count_elem, &h_);

	int i_count[6];
	i_count[0] = 0;
	double x[6] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };
	int t = count_elem / 5;
	for (int j = 1; j < 6; j++) {
		i_count[j] = t * j;
	}
	i_count[5]--;

	double ta = a;
	for (int i = 0; i < count_elem; i++, ta += h_) {
		fprintf(out, "%.4lf %lf\n", ta, yt[i][0]);
	}

	vector<pair<double, double>> v;

	double m = 0.0;
	for (int i = 0; i < 6; i++, m += h) {
		printf("%.2lf\t", m);
		printf("%.3lf\t", yt[i_count[i]][0]);
		printf("%.3lf\n", yt[i_count[i]][1]);
		y[i] = yt[i_count[i]][0];

		v.push_back(make_pair(x[i], y[i]));
	}
	printf("\n");

	fclose(out);

	printf("Spleins interpolation in file\n");
	FILE *splines_out = fopen("Splines.txt", "w");
	

	for (double i = a; i <= b; i += h_) {
		double tmp = cubic(i, v.size(), v.data());
		fprintf(splines_out, "%lf %lf\n", i, tmp);
	}
	printf("\n");

	fclose(splines_out);

	Eps = 1e-2;
	// printf("Integration = %.10lf\n", Form_of_Simpson(a, b, h, tmp));
	printf("Integration = %.10lf\n", double_counting(Form_of_Simpson, a, b, h_, Eps, tmp));

	return 0;
}