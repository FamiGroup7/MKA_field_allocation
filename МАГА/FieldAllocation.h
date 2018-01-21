#pragma once
#include "Mka2D_cylindrical.h"
#include "Mka3D.h"

#define NOPE false

using namespace std;

class FIeldAllocation {
public:
	FIeldAllocation(double lambda0, double power, bool withCheck);
	~FIeldAllocation();

	void start();
	void start_without2D();

	void output_slice(const string & tecplot_filename, char slice_var, double slice_val, char var1, double min_var1, double max_var1, size_t num_var_1, char var2, double min_var2, double max_var2, size_t num_var_2);

	void output_slice(const string & tecplot_filename, double *q0, char slice_var, double slice_val, char var1, double min_var1, double max_var1, size_t num_var_1, char var2, double min_var2, double max_var2, size_t num_var_2);
	void output_slice_sameNet(const string & tecplot_filename, double * q0, char slice_var, double slice_val, char var1, double min_var1, double max_var1, size_t num_var_1, char var2, double min_var2, double max_var2, size_t num_var_2);

	double *q0;
private:
	bool withCheck;
	Point source;
	double power;
	double lambda0;
	Mka2D_cylindrical*task2D;
	Mka3D*taskSimple;
	Mka3D*taskCheck;
	ofstream logger;

	double * map2dSolutionToNet(Mka2D_cylindrical & from, Mka3D & to);
	void generageG_lambda(double lambda0, Mka3D & simple);
	void fieldAllocation(double lambda0, double * q0, Mka3D & simple);
	double * analitic_u0(double power, Point source, double lambda, Mka3D & net);

	struct data_elem
	{
		double v1, v2;
		double u_correct, u0, u_plus;
	};
};