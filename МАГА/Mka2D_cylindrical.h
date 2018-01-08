#pragma once
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <Windows.h>
#include <cstring>
#include <math.h>
#include <algorithm>
#include "Point_cylindrical.h"
#include "Point.h"

using namespace std;

class Mka2D_cylindrical
{
public:
	Mka2D_cylindrical();
	Mka2D_cylindrical(string filePrefix, bool ku1_left, bool ku1_right, bool ku1_up, bool ku1_down);
	void startFull2dProcess();
	~Mka2D_cylindrical();

	int FindAreaNumber(int nodes[]);
	void GenerateNet();
	void GeneratePortrait();
	double analiticSolution(Point_cylindrical source);
	double Func(int ielem, int i);
	double BasicFunc1d(int num, double ksi);
	double BasicFunc2d(int i, double ksi, double eta);
	double difBasicFunc1d(int num, double ksi);
	void CreateLocalMatrixs(int ielem);
	void Addition(int ielem);
	void PrintPlotMatrix(bool flag_simmeric);
	int indexRZ(Point_cylindrical source);
	int indexRZ(double r, double z);
	void Edge1_sim();
	void Edge1_sim_old();
	void Edge2();
	void GenerateMatrix();
	void genNet1d(double startValue, double endValue, double startH, double koefRazriadki, vector<double>& vect);
	void insertInVector(double value, vector<double>& vect);
	void MultMatrixOnVector(double * in, double * out);
	void MultMatrixOnVector(double * in, double * out, double * diMas, double * gglMas, double * gguMas);
	double ScalarMult(double * v1, double * v2);
	void runLOS(double * ggl, double * ggu, double * diag, int N, int * ig, int * jg, double * b, double * q);
	void LOS();
	void directSolveStraightTask();
	int findKE(Point_cylindrical point);

	double power = 0;

	double solutionInPoint(int iKe, Point_cylindrical target);

	struct field
	{
		double r1, r2, z1, z2, sigma, gamma;
		field(double x1_new, double x2_new, double y1_new, double y2_new, double lambda_new, double gamma_new)
		{
			r1 = x1_new;
			r2 = x2_new;
			z1 = y1_new;
			z2 = y2_new;
			sigma = lambda_new;
			gamma = gamma_new;
		}
		field()
		{
			r1 = r2 = z1 = z2 = sigma = gamma = 0;
		}
	};

	struct nvtr
	{
		nvtr(int uz1, int uz2, int uz3, int uz4, int numberField_new)
		{
			uzel[0] = uz1; uzel[1] = uz2; uzel[2] = uz3; uzel[3] = uz4;
			numberField = numberField_new;
		}
		nvtr(int uz1, int uz2, int uz3, int uz4)
		{
			uzel[0] = uz1; uzel[1] = uz2; uzel[2] = uz3; uzel[3] = uz4;
			numberField = 0;
		}
		nvtr()
		{
			uzel[0] = uzel[1] = uzel[2] = uzel[3] = numberField = 0;
		}
		int uzel[4], numberField;
	};

	double* q;
	vector<double> R;
	vector<double> Z;
	vector<field> sreda;
	vector<Point_cylindrical> rz;
	vector<nvtr> KE;
private:
	string filePrefix;
	ofstream testFile;
	int nKE, nPoints;
	int *ig, *jg;
	double *ggl, *ggu, *di, *b;

	double e = 1e-16;
	int maxiter = 10000;

	double localB[4], localMatrix[4][4];
	double helpG1[4][4] = { { 2, -2, 1, -1 },{ -2, 2, -1, 1 },{ 1, -1, 2, -2 },{ -1, 1, -2, 2 } };
	double helpG2[4][4] = { { 2, 1, -2, -1 },{ 1, 2, -1, -2 },{ -2, -1, 2, 1 },{ -1, -2, 1, 2 } };
	double helpGOseSim[4][4] = { { 1, 1, -1, -1 },{ 1, 3, -1, -3 },{ -1, -1, 1, 1 },{ -1, -3, 1, 3 } };
	double helpM[4][4] = { { 4, 2, 2, 1 },{ 2, 4, 1, 2 },{ 2, 1, 4, 2 },{ 1, 2, 2, 4 } };

	double r0, r1, z0, z1; 
	bool ku1_left, ku1_right, ku1_up, ku1_down;
	double koordSourceR, koordSourceZ;
};

