#include <math.h>
#include <iostream>

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
