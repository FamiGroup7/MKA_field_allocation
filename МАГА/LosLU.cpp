#pragma once
#include <math.h>
#include "MkaUtils.h"

using namespace std;

void LUsq(double* ggl, double* ggu, double* diag, int N, int* ig, int* jg, double* D, double* L, double* U)
{
	double sd;
	for (int i = 0; i < N; i++)
	{
		D[i] = diag[i];
		sd = 0;
		for (int j = ig[i], end = ig[i + 1]; j < end; j++)
		{
			double sl = 0, su = 0;
			int k = jg[j];
			int on_k = ig[k];
			int on_i = ig[i];
			int end_k = ig[k + 1];

			while (on_k < end_k && on_i < j)
			{
				if (jg[on_i] == jg[on_k])
				{
					sl += L[on_i] * U[on_k];
					su += U[on_i] * L[on_k];
					on_i++;
					on_k++;
				}
				else if (jg[on_i] < jg[on_k]) on_i++;
				else on_k++;
			}

			L[j] = (ggl[j] - sl) / D[k];
			U[j] = (ggu[j] - su) / D[k];
			sd += U[j] * L[j];
		}
		if(D[i]<sd)
		{
			std::cerr << "Error in LUsq. D[i]<sd" << std::endl;
			system("pause");
			exit(1);
		}
		D[i] = sqrt(D[i] - sd);
	}
}

void mult(double* a1, double* a2, double* ggl, double* ggu, double* diag, int N, int* ig, int* jg)
{
	for (int i = 0; i < N; i++)
	{
		a1[i] = 0;
	}

	for (int i = 0; i < N; i++)
	{
		a1[i] += diag[i] * a2[i];
		for (int j = ig[i], end = ig[i + 1]; j < end; j++)
		{
			a1[i] += ggl[j] * a2[jg[j]];
			a1[jg[j]] += ggu[j] * a2[i];
		}
	}
}

double scalarOfVectors(double* a1, double* a2, int dim)
{
	double S = 0.0;

	for (int i = 0; i < dim; i++)
		S += a1[i] * a2[i];

	return S;
}

void multL(double* a1, int N, int* ig, int* jg, double* L, double* D)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = ig[i], end = ig[i + 1]; j < end; j++)
		{
			a1[i] -= L[j] * a1[jg[j]];
		}

		a1[i] = a1[i] / D[i];
	}
}

void multU(double* a1, int N, int* ig, int* jg, double* U, double* D)
{
	for (int i = N - 1; i >= 0; i--)
	{
		a1[i] = a1[i] / D[i];
		for (int j = ig[i], end = ig[i + 1]; j < end; j++)
		{
			a1[jg[j]] -= U[j] * a1[i];
		}
	}
}

int LOS(int n, int *ig, int *jg,
	double *ggl, double *ggu, double *di, double *f,
	double *x, double eps, int maxiter)
{
	int iter = 0;
	int i;
	double g, g1, mu, alpha, betta;
	double nu;
	double *r = NULL, *p = NULL, *h = NULL, *s = NULL, *z = NULL, *tmp = NULL, *w = NULL;
	double r_old, p_scal;

	__time64_t time_total, time_beg, time_end; // для засечки времени
	__time64_t timestruct;

	_time64(&timestruct);
	time_beg = timestruct;

	ofstream logfile("solver.txt");
	logfile.open("log.txt");
	logfile << "LOS..." << endl;

	//// факторизованнные диагональные блоки
	//double *df = NULL;     // диагональ факторизованных диагональных блоков 
	//double *ggl_f = NULL;  // внедиагональные эл-ты факторизованных диагональных блоков (нижний треугольник)
	//double *ggu_f = NULL;  // внедиагональные эл-ты факторизованных диагональных блоков (верхний треугольник)

	r = new double[n];
	p = new double[n];
	s = new double[n];
	z = new double[n];
	tmp = new double[n];
	w = new double[n];
	h = new double[n];

	//double *s_smooth = new double[n];
	//double *y = new double[n];

	//df = new double[n];
	//ggl_f = new double[n];
	//ggu_f = new double[n];

	//// выполняем факторизацию диагональных блоков
	//Build_block_diag_preconditioner(nb, idi, di_block, df, ggl_f, ggu_f);

	// начальное приближение
	for (i = 0; i < n; i++)
	{
		r[i] = h[i] = s[i] = z[i] = x[i] = p[i] = 0.0;
	}

	//	вычисление начальной невязки
	mult(r, x, ggl, ggu, di, n, ig, jg);

	for (i = 0; i < n; i++)
	{
		r[i] = f[i] - r[i];
		p[i] = r[i];
	}
	mult(p, z, ggl, ggu, di, n, ig, jg);

	betta = 0;

	r_old = scalarOfVectors(r, r, n);

	//for (size_t i = 0; i < n; i++)
	//{
	//	s_smooth[i] = r[i];
	//	y[i] = x[i];
	//	p[i] = s[i];
	//}
	//s_old = scalarOfVectors(s_smooth, s_smooth, n);

	if (r_old < 1e-30)
	{
		logfile << "x0 is solution" << endl;
		//cout << "x0 is solution" << endl;
		goto COCG_EXIT;
	}

	// главный цикл
	for (iter = 1; iter <= maxiter; iter++)
	{
		g = scalarOfVectors(p, r, n);
		p_scal = scalarOfVectors(p, p, n);
		alpha = g / p_scal;

		// x = x + alpha*p
		for (i = 0; i < n; i++)
			x[i] += alpha*z[i];

		// r = r - alpha*u
		for (i = 0; i < n; i++)
			r[i] -= alpha*p[i];

		mult(tmp, r, ggl, ggu, di, n, ig, jg);

		g1 = scalarOfVectors(p, tmp, n);
		betta = -g1 / p_scal;

		// p = s + betta*p
		for (i = 0; i < n; i++)
			z[i] = r[i] + betta*z[i];

		// z = tmp + betta*z
		for (i = 0; i < n; i++)
			p[i] = tmp[i] + betta*p[i];

		if (fabs(p_scal) < 1e-30 && fabs(p_scal) < 1e-30)
		{
			//cout << "mu=0, LOS failed" << endl;
			logfile << "mu=0, LOS failed" << endl;
			goto COCG_EXIT;
		}


		//// nu

		//double *rs = new double[n];
		//for (size_t i = 0; i < n; i++)
		//{
		//	rs[i] = r[i] - s_smooth[i];
		//}

		//nu = Scal(s_smooth, rs, n) / Scal(rs, rs, n) * (-1.0);
		//if (nu < 0.0)
		//	nu = 0.0;
		//if (nu > 1.0)
		//	nu = 1.0;

		//for (i = 0; i < n; i++)
		//{
		//	y[i] = (1.0 - nu) * y[i] + nu * x[i];
		//	s_smooth[i] = (1.0 - nu) * s_smooth[i] + nu * r[i];
		//}


		// вычисляем норму невязки
		//r = Norm_Euclid(r, n);
		double s_norm = scalarOfVectors(r, r, n);

		//cout << iter << "\t" << scientific << r / p_old << endl;
		logfile << iter << "\t" << scientific << s_norm / r_old << endl;
		//logfile << iter << "\t" << scientific << r / p_old << endl;


		// выход по достижении малости невязки
		//if (r / p_old < eps)
		//	goto COCG_EXIT;
		if (s_norm / r_old < eps)
			goto COCG_EXIT;
		//p_old = r;


		if (fabs(g1) < 1e-30 && fabs(g1) < 1e-30)
		{
			//cout << "g=0, LOS failed" << endl;
			logfile << "g=0, LOS failed" << endl;
			goto COCG_EXIT;
		}
		//if (rs) { delete[] rs; rs = NULL; }
	}

COCG_EXIT:

	cout << "iter: " << iter << endl;
	if (r) { delete[] r; r = NULL; }
	if (p) { delete[] p; p = NULL; }
	//if (s_smooth) { delete[] s_smooth; s_smooth = NULL; }
	if (h) { delete[] h; h = NULL; }
	if (s) { delete[] s; s = NULL; }
	//if (y) { delete[] y; y = NULL; }

	//if (df) { delete[] df; df = NULL; }
	//if (ggl_f) { delete[] ggl_f; ggl_f = NULL; }
	//if (ggu_f) { delete[] ggu_f; ggu_f = NULL; }

	_time64(&timestruct);
	time_end = timestruct;
	time_total = time_end - time_beg;
	logfile << "Time: " << time_total << endl;
	//	Write_kit("kit", r / p_old, eps, iter, time_total);
	logfile.close();

	return iter;
}

int LosLU(double* ggl, double* ggu, double* diag, int N, int* ig, int* jg, double* f, double* q)
{
	double normR = 0.0, normR0 = 0.0, alpha = 0.0, beta = 0.0;
	int maxiter = 10000;
	double eps = 1e-16;
	//clear(A.N, A.ig);
	double* D = new double[N];
	double* L = new double[ig[N]];
	double* U = new double[ig[N]];
	double* ri = new double[N];
	double* zi = new double[N];
	double* pi = new double[N];
	double* Utek = new double[N];
	double* tek = new double[N];
	LUsq(ggl, ggu, diag, N, ig, jg, D, L, U);
	for (int i = 0; i < N; i++)
	{
		q[i] = zi[i] = pi[i] = Utek[i] = tek[i] = 0;
		ri[i] =  f[i];
	}
	//	нашли r0
	multL(ri, N, ig, jg, L, D);
	//	посчитали норму невязки
	normR0 = scalarOfVectors(ri, ri, N);
	normR0 = sqrt(normR0);
	//	установили начальные значения векторов
	//	нашли z0
	for (int i = 0; i < N; i++)
	{
		zi[i] = ri[i];
	}
	multU(zi, N, ig, jg, U, D);
	//	нашли p0
	mult(pi, zi, ggl, ggu, diag, N, ig, jg);
	multL(pi, N, ig, jg, L, D);
	double nev = 0.0;

	int k = 1;
	do
	{
		//	находим альфу
		alpha = scalarOfVectors(pi, ri, N) / scalarOfVectors(pi, pi, N);
		normR = scalarOfVectors(ri, ri, N) - alpha * alpha * scalarOfVectors(pi, pi, N);
		normR = sqrt(normR);
		nev = normR / normR0;

		for (int i = 0; i < N; ++i)
		{
			q[i] += alpha * zi[i];
			ri[i] -= alpha * pi[i];
		}

		//	находим вектор, чтобы пересчитать z, p
		for (int i = 0; i < N; i++)
		{
			Utek[i] = ri[i];
		}
		multU(Utek, N, ig, jg, U, D);
		mult(tek, Utek, ggl, ggu, diag, N, ig, jg);
		multL(tek, N, ig, jg, L, D);

		beta = - scalarOfVectors(pi, tek, N) / scalarOfVectors(pi, pi, N);

		for (int i = 0; i < N; ++i)
		{
			zi[i] = Utek[i] + beta * zi[i];
			pi[i] = tek[i] + beta * pi[i];
		}
		k++;
	}
	while (k < maxiter && nev >= eps);

	//printf("	Iter = %d; nev = %le\n", k, nev);
	return k;
}
