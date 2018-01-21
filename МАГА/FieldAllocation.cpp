#include "FieldAllocation.h"

int LosLU(double* ggl, double* ggu, double* diag, int N, int* ig, int* jg, double* f, double* q);

FIeldAllocation::FIeldAllocation(double lambda0, double power, bool withCheck)
	: lambda0(lambda0), power(power), withCheck(withCheck)
{
	logger.open("field_allocation.log");
	if (withCheck == true) {
		taskCheck = new Mka3D("resources_field_allocation/3D_check/", NOPE, NOPE, NOPE, NOPE, NOPE, NOPE, NOPE);
	}
	task2D = new Mka2D_cylindrical("resources_field_allocation/2D/", false, true, false, true);
	taskSimple = new Mka3D("resources_field_allocation/3D/", NOPE, NOPE, NOPE, NOPE, NOPE, NOPE, NOPE);

	source.x = 0; source.y = 0; source.z = 0;
	task2D->koordSourceR = 0;
	task2D->koordSourceZ = 0;
	taskSimple->koordSourceX = 0;
	taskSimple->koordSourceY = 0;
	taskSimple->koordSourceZ = 0;
	taskCheck->koordSourceX = 0;
	taskCheck->koordSourceY = 0;
	taskCheck->koordSourceZ = 0;
}

FIeldAllocation::~FIeldAllocation()
{
	delete task2D;
	delete taskSimple;
	delete taskCheck;
}

void FIeldAllocation::start()
{
	logger << "Field allocation" << endl;
	auto start = std::chrono::system_clock::now();
	task2D->GenerateNet();
	task2D->power = power;
	task2D->sourceType = 1;
	task2D->directSolveStraightTask();
	logger << "2D task" << endl;
	logger << "KE.size=" << task2D->KE.size() << ", rz.size=" << task2D->rz.size() << endl;
	logger << setw(40) << std::left << "duration" <<
		MkaUtils::formattingTime(std::chrono::system_clock::time_point(std::chrono::system_clock::now() - start)) << endl;

	start = std::chrono::system_clock::now();
	taskSimple->buildNet("sreda.txt", "absent");
	taskSimple->build_xyz_nvtr_portratin_Abqx();
	taskSimple->netOptimization();
	taskSimple->power = power;
	taskSimple->sourceType = 1;
	logger << "construct simple 3D net" << endl;
	logger << "KE.size=" << taskSimple->KE.size() << ", xyz.size=" << taskSimple->xyz_points.size() << endl;
	logger << setw(40) << std::left << "duration" <<
		MkaUtils::formattingTime(std::chrono::system_clock::time_point(std::chrono::system_clock::now() - start)) << endl;

	start = std::chrono::system_clock::now();
	q0 = map2dSolutionToNet(*task2D, *taskSimple);
	logger << setw(40) << std::left << "map 2D to 3D simple net" <<
		MkaUtils::formattingTime(std::chrono::system_clock::time_point(std::chrono::system_clock::now() - start)) << endl;

	start = std::chrono::system_clock::now();
	fieldAllocation(lambda0, q0, *taskSimple);
	logger << setw(40) << std::left << "calc q+" <<
		MkaUtils::formattingTime(std::chrono::system_clock::time_point(std::chrono::system_clock::now() - start)) << endl;

	//check
	if (!withCheck)return;
	start = std::chrono::system_clock::now();
	taskCheck->buildNet("sreda.txt", "sourceLocate.txt");
	taskCheck->netOptimization();
	taskCheck->build_xyz_nvtr_portratin_Abqx();
	taskCheck->power = power;
	taskCheck->sourceType = 1;
	taskCheck->generateGlobalMatrix();
	LosLU(taskCheck->ggl, taskCheck->ggu, taskCheck->di, taskCheck->countRegularNodes, 
		taskCheck->ig, taskCheck->jg, taskCheck->b, taskCheck->q);
	logger << "KE.size=" << taskCheck->KE.size() << ", xyz.size=" << taskCheck->xyz_points.size() << endl;
	logger << setw(40) << std::left << "check task duration" <<
		MkaUtils::formattingTime(std::chrono::system_clock::time_point(std::chrono::system_clock::now() - start)) << endl;
}

void FIeldAllocation::start_without2D() {
	logger << "Field allocation without 2D" << endl;

	auto start = std::chrono::system_clock::now();
	taskSimple->buildNet("sreda.txt", "absent");
	taskSimple->build_xyz_nvtr_portratin_Abqx();
	taskSimple->netOptimization();
	taskSimple->power = power;
	taskSimple->sourceType = 1;
	taskSimple->generateGlobalMatrix(lambda0);
	LosLU(taskSimple->ggl, taskSimple->ggu, taskSimple->di, taskSimple->countRegularNodes,
		taskSimple->ig, taskSimple->jg, taskSimple->b, taskSimple->q);
	logger << "solve simple 3D with only lambda0" << endl;
	logger << "KE.size=" << taskSimple->KE.size() << ", xyz.size=" << taskSimple->xyz_points.size() << endl;
	logger << setw(40) << std::left << "duration" <<
		MkaUtils::formattingTime(std::chrono::system_clock::time_point(std::chrono::system_clock::now() - start)) << endl;

	start = std::chrono::system_clock::now();
	q0 = new double[taskSimple->countRegularNodes];
	memcpy(q0, taskSimple->q, sizeof(double)*taskSimple->countRegularNodes);
	fieldAllocation(lambda0, q0, *taskSimple);
	logger << setw(40) << std::left << "calc q+" <<
		MkaUtils::formattingTime(std::chrono::system_clock::time_point(std::chrono::system_clock::now() - start)) << endl;

	//check
	if (!withCheck)return;
	start = std::chrono::system_clock::now();
	taskCheck->buildNet("sreda.txt", "sourceLocate.txt");
	taskCheck->netOptimization();
	taskCheck->build_xyz_nvtr_portratin_Abqx();
	taskCheck->power = power;
	taskCheck->sourceType = 1;
	taskCheck->generateGlobalMatrix();
	LosLU(taskCheck->ggl, taskCheck->ggu, taskCheck->di, taskCheck->countRegularNodes,
		taskCheck->ig, taskCheck->jg, taskCheck->b, taskCheck->q);
	logger << "KE.size=" << taskCheck->KE.size() << ", xyz.size=" << taskCheck->xyz_points.size() << endl;
	logger << setw(40) << std::left << "check task duration" <<
		MkaUtils::formattingTime(std::chrono::system_clock::time_point(std::chrono::system_clock::now() - start)) << endl;
}

double* FIeldAllocation::map2dSolutionToNet(Mka2D_cylindrical&from, Mka3D&to) {
	double*q0 = new double[to.countRegularNodes];
	for (int i = 0; i < to.countRegularNodes; i++)
	{
		Point target = to.xyz_points[i];
		Point_cylindrical target2D(sqrt(target.x*target.x + target.y*target.y), target.z);
		int iKe2D = from.findKE(target2D);
		double solution = 0;
		if (iKe2D >= 0) {
			solution = from.solutionInPoint(iKe2D, target2D);
		}
		q0[i] = solution;
	}
	return q0;
}

void FIeldAllocation::generageG_lambda(double lambda0, Mka3D&simple) {
	memset(simple.di, 0, sizeof(double)*simple.countRegularNodes);
	memset(simple.ggl, 0, sizeof(double)*simple.ig[simple.countRegularNodes]);

	//create G_delta_lambda
	for (int ielem = 0; ielem < simple.KE.size(); ielem++)
	{
		int area = simple.FindAreaNumber(simple.KE[ielem].uzel);
		double val = lambda0 - simple.sreda[area].lambda;
		if (MkaUtils::equals(val, 0)) {
			continue;
		}
		simple.CreateLocalMatrixs(ielem, val);
		simple.Addition(ielem, simple.di, simple.ggl);
	}

	for (int i = 0; i < simple.ig[simple.countRegularNodes]; i++)
	{
		simple.ggu[i] = simple.ggl[i];
	}
}

//return solution u0
void FIeldAllocation::fieldAllocation(double lambda0, double*q0, Mka3D&simple) {
	generageG_lambda(lambda0, simple);

	double* b_plus = new double[simple.countRegularNodes];
	memset(b_plus, 0, sizeof(double)*simple.countRegularNodes);
	simple.MultMatrixOnVector(q0, b_plus, simple.di, simple.ggl, simple.ggu);

	for (int i = 0; i < simple.ig[simple.countRegularNodes]; i++)
	{
		simple.ggu[i] = simple.ggl[i] = 0;
	}
	for (int i = 0; i < simple.countRegularNodes; i++) {
		simple.di[i] = simple.b[i] = simple.q[i] = 0;
	}
	//1e kraevoe
	simple.leftKU = simple.rightKU/* = simple.upKU*/ = simple.downKU = simple.foreKU = simple.behindKU = 1;
	simple.upKU = 5;
	simple.power = 0;
	simple.sourceType = 0;
	simple.generateGlobalMatrix(/*lambda0*/);
	LosLU(simple.ggl, simple.ggu, simple.di, simple.countRegularNodes, simple.ig, simple.jg, b_plus, simple.q);
}

double* FIeldAllocation::analitic_u0(double power, Point source, double lambda, Mka3D&net) {
	double koef = power / (2 * PI/**lambda*/);
	double *q0 = new double[net.countRegularNodes];
	for (size_t i = 0; i < net.countRegularNodes; i++)
	{
		q0[i] = koef / (sqrt(
			(source.x - net.xyz_points[i].x)*(source.x - net.xyz_points[i].x) +
			(source.y - net.xyz_points[i].y)*(source.y - net.xyz_points[i].y) +
			(source.z - net.xyz_points[i].z)*(source.z - net.xyz_points[i].z)
		));
	}
	return q0;
}

//--------------

void FIeldAllocation::output_slice(const string & tecplot_filename,
	char slice_var, double slice_val,
	char var1, double min_var1, double max_var1, size_t num_var_1,
	char var2, double min_var2, double max_var2, size_t num_var_2)
{
	if (min_var1 > max_var1) swap(max_var1, min_var1);
	if (min_var2 > max_var2) swap(max_var2, min_var2);

	double step_var1 = (max_var1 - min_var1) / (double)(num_var_1 - 1);
	double step_var2 = (max_var2 - min_var2) / (double)(num_var_2 - 1);
	if (num_var_1 < 2)step_var1 = 0;
	if (num_var_2 < 2)step_var2 = 0;
	size_t index_slice = 0, index_1 = 0, index_2 = 0;

	// Определяем по какой переменной сечение
	if (slice_var == 'x' || slice_var == 'X') index_slice = 0;
	else if (slice_var == 'y' || slice_var == 'Y') index_slice = 1;
	else if (slice_var == 'z' || slice_var == 'Z') index_slice = 2;
	else
	{
		cerr << "Unknown variable, breaking ..." << endl;
		return;
	}

	// Определяем, какая переменная первая
	if (var1 == 'x' || var1 == 'X') index_1 = 0;
	else if (var1 == 'y' || var1 == 'Y') index_1 = 1;
	else if (var1 == 'z' || var1 == 'Z') index_1 = 2;
	else
	{
		cerr << "Unknown variable, breaking ..." << endl;
		return;
	}

	// Определяем, какая переменная вторая
	if (var2 == 'x' || var2 == 'X') index_2 = 0;
	else if (var2 == 'y' || var2 == 'Y') index_2 = 1;
	else if (var2 == 'z' || var2 == 'Z') index_2 = 2;
	else
	{
		cerr << "Unknown variable, breaking ..." << endl;
		return;
	}

	cout << "Writing slice to Tecplot ..." << endl;

	ofstream tecplot_file;
	tecplot_file.open(tecplot_filename.c_str(), ios::out);

	if (!tecplot_file.good())
	{
		cerr << "Error in " << __FILE__ << ":" << __LINE__
			<< " while writing file " << tecplot_filename << endl;
		return;
	}

	tecplot_file << "TITLE = \"Slice " << slice_var << " = " << slice_val << "\"\n";
	tecplot_file << "VARIABLES = \"" << var1 << "\", \"" << var2 << "\", \"V\", \"Ex\", \"Ey\", \"Ez\", \"abs(E)\"\n";
	tecplot_file << "ZONE I= " << num_var_1 << ", J= " << num_var_2 << ", F=POINT\n";

	tecplot_file.precision(17);
	tecplot_file.setf(ios::scientific);

	vector<data_elem> data(num_var_1 * num_var_2);

	for (size_t i = 0; i < num_var_1; i++)
	{
		Point p(0, 0, 0);
		p[index_slice] = slice_val;
		data_elem tmp_elem;
		double v1 = min_var1 + step_var1 * (double)i;
		tmp_elem.v1 = v1;
		for (size_t j = 0; j < num_var_2; j++)
		{
			double v2 = min_var2 + step_var2 * (double)j;
			tmp_elem.v2 = v2;
			p[index_1] = v1;
			p[index_2] = v2;
			int iKE_allocation = taskSimple->findKE(p);
			if (iKE_allocation < 0) {
				throw new exception();
			}
			tmp_elem.u_plus = taskSimple->solutionInPoint(iKE_allocation, p);
			tmp_elem.u0 = taskSimple->solutionInPoint(iKE_allocation, p, q0);

			int iKE = taskCheck->findKE(p);
			if (iKE < 0) {
				throw new exception();
			}
			tmp_elem.u_correct = taskCheck->solutionInPoint(iKE, p);

			data[i * num_var_2 + j] = tmp_elem;
		}
	}

	tecplot_file << "v1;v2;u0;u_plus;u_correct;u_correct-(u0+u_plus)" << endl;
	for (size_t i = 0; i < data.size(); i++)
		tecplot_file << data[i].v1 << ';' << data[i].v2 << ';' << data[i].u0 << ';' << data[i].u_plus
		<< ';' << data[i].u_correct << ';' << data[i].u_correct - (data[i].u0 + data[i].u_plus) << endl;

	tecplot_file << "\n";
	tecplot_file.flush();
	tecplot_file.close();
}

void FIeldAllocation::output_slice(const string & tecplot_filename, double*q0,
	char slice_var, double slice_val,
	char var1, double min_var1, double max_var1, size_t num_var_1,
	char var2, double min_var2, double max_var2, size_t num_var_2)
{
	if (min_var1 > max_var1) swap(max_var1, min_var1);
	if (min_var2 > max_var2) swap(max_var2, min_var2);

	double step_var1 = (max_var1 - min_var1) / (double)(num_var_1 - 1);
	double step_var2 = (max_var2 - min_var2) / (double)(num_var_2 - 1);
	if (num_var_1 < 2)step_var1 = 0;
	if (num_var_2 < 2)step_var2 = 0;
	size_t index_slice = 0, index_1 = 0, index_2 = 0;

	// Определяем по какой переменной сечение
	if (slice_var == 'x' || slice_var == 'X') index_slice = 0;
	else if (slice_var == 'y' || slice_var == 'Y') index_slice = 1;
	else if (slice_var == 'z' || slice_var == 'Z') index_slice = 2;
	else
	{
		cerr << "Unknown variable, breaking ..." << endl;
		return;
	}

	// Определяем, какая переменная первая
	if (var1 == 'x' || var1 == 'X') index_1 = 0;
	else if (var1 == 'y' || var1 == 'Y') index_1 = 1;
	else if (var1 == 'z' || var1 == 'Z') index_1 = 2;
	else
	{
		cerr << "Unknown variable, breaking ..." << endl;
		return;
	}

	// Определяем, какая переменная вторая
	if (var2 == 'x' || var2 == 'X') index_2 = 0;
	else if (var2 == 'y' || var2 == 'Y') index_2 = 1;
	else if (var2 == 'z' || var2 == 'Z') index_2 = 2;
	else
	{
		cerr << "Unknown variable, breaking ..." << endl;
		return;
	}

	cout << "Writing slice to Tecplot ..." << endl;

	ofstream tecplot_file;
	tecplot_file.open(tecplot_filename.c_str(), ios::out);

	if (!tecplot_file.good())
	{
		cerr << "Error in " << __FILE__ << ":" << __LINE__
			<< " while writing file " << tecplot_filename << endl;
		return;
	}

	tecplot_file << "TITLE = \"Slice " << slice_var << " = " << slice_val << "\"\n";
	tecplot_file << "VARIABLES = \"" << var1 << "\", \"" << var2 << "\", \"V\", \"Ex\", \"Ey\", \"Ez\", \"abs(E)\"\n";
	tecplot_file << "ZONE I= " << num_var_1 << ", J= " << num_var_2 << ", F=POINT\n";

	tecplot_file.precision(17);
	tecplot_file.setf(ios::scientific);

	vector<data_elem> data(num_var_1 * num_var_2);

	for (size_t i = 0; i < num_var_1; i++)
	{
		Point p(0, 0, 0);
		p[index_slice] = slice_val;
		data_elem tmp_elem;
		double v1 = min_var1 + step_var1 * (double)i;
		tmp_elem.v1 = v1;
		for (size_t j = 0; j < num_var_2; j++)
		{
			double v2 = min_var2 + step_var2 * (double)j;
			tmp_elem.v2 = v2;
			p[index_1] = v1;
			p[index_2] = v2;
			int iKE_allocation = taskSimple->findKE(p);
			if (iKE_allocation < 0) {
				throw new exception();
			}
			tmp_elem.u_plus = taskSimple->solutionInPoint(iKE_allocation, p);
			tmp_elem.u0 = taskSimple->solutionInPoint(iKE_allocation, p, q0);

			int iKE = taskCheck->findKE(p);
			if (iKE < 0) {
				throw new exception();
			}
			tmp_elem.u_correct = taskCheck->solutionInPoint(iKE, p);

			data[i * num_var_2 + j] = tmp_elem;
		}
	}

	tecplot_file << "v1;v2;u0;u_plus;u_correct;u_correct-(u0+u_plus)" << endl;
	for (size_t i = 0; i < data.size(); i++)
		tecplot_file << data[i].v1 << ';' << data[i].v2 << ';' << data[i].u0 << ';' << data[i].u_plus
		<< ';' << data[i].u_correct << ';' << data[i].u_correct - (data[i].u0 + data[i].u_plus) << endl;

	tecplot_file << "\n";
	tecplot_file.flush();
	tecplot_file.close();
}

//--------------------------------------

void FIeldAllocation::output_slice_sameNet(const string & tecplot_filename, double*q0,
	char slice_var, double slice_val,
	char var1, double min_var1, double max_var1, size_t num_var_1,
	char var2, double min_var2, double max_var2, size_t num_var_2)
{
	if (min_var1 > max_var1) swap(max_var1, min_var1);
	if (min_var2 > max_var2) swap(max_var2, min_var2);

	double step_var1 = (max_var1 - min_var1) / (double)(num_var_1 - 1);
	double step_var2 = (max_var2 - min_var2) / (double)(num_var_2 - 1);
	if (num_var_1 < 2)step_var1 = 0;
	if (num_var_2 < 2)step_var2 = 0;
	size_t index_slice = 0, index_1 = 0, index_2 = 0;

	// Определяем по какой переменной сечение
	if (slice_var == 'x' || slice_var == 'X') index_slice = 0;
	else if (slice_var == 'y' || slice_var == 'Y') index_slice = 1;
	else if (slice_var == 'z' || slice_var == 'Z') index_slice = 2;
	else
	{
		cerr << "Unknown variable, breaking ..." << endl;
		return;
	}

	// Определяем, какая переменная первая
	if (var1 == 'x' || var1 == 'X') index_1 = 0;
	else if (var1 == 'y' || var1 == 'Y') index_1 = 1;
	else if (var1 == 'z' || var1 == 'Z') index_1 = 2;
	else
	{
		cerr << "Unknown variable, breaking ..." << endl;
		return;
	}

	// Определяем, какая переменная вторая
	if (var2 == 'x' || var2 == 'X') index_2 = 0;
	else if (var2 == 'y' || var2 == 'Y') index_2 = 1;
	else if (var2 == 'z' || var2 == 'Z') index_2 = 2;
	else
	{
		cerr << "Unknown variable, breaking ..." << endl;
		return;
	}

	cout << "Writing slice to Tecplot ..." << endl;

	ofstream tecplot_file;
	tecplot_file.open(tecplot_filename.c_str(), ios::out);

	if (!tecplot_file.good())
	{
		cerr << "Error in " << __FILE__ << ":" << __LINE__
			<< " while writing file " << tecplot_filename << endl;
		return;
	}

	tecplot_file << "TITLE = \"Slice " << slice_var << " = " << slice_val << "\"\n";
	tecplot_file << "VARIABLES = \"" << var1 << "\", \"" << var2 << "\", \"V\", \"Ex\", \"Ey\", \"Ez\", \"abs(E)\"\n";
	tecplot_file << "ZONE I= " << num_var_1 << ", J= " << num_var_2 << ", F=POINT\n";

	tecplot_file.precision(17);
	tecplot_file.setf(ios::scientific);

	vector<data_elem> data(num_var_1 * num_var_2);

	for (size_t i = 0; i < num_var_1; i++)
	{
		Point p(0, 0, 0);
		p[index_slice] = slice_val;
		data_elem tmp_elem;
		double v1 = min_var1 + step_var1 * (double)i;
		tmp_elem.v1 = v1;
		for (size_t j = 0; j < num_var_2; j++)
		{
			double v2 = min_var2 + step_var2 * (double)j;
			tmp_elem.v2 = v2;
			p[index_1] = v1;
			p[index_2] = v2;
			int iKE = taskSimple->findKE(p);
			if (iKE < 0) {
				throw new exception();
			}
			tmp_elem.u_correct = taskSimple->solutionInPoint(iKE, p);
			tmp_elem.u0 = taskSimple->solutionInPoint(iKE, p, q0);
			data[i * num_var_2 + j] = tmp_elem;
		}
	}

	tecplot_file << "v1;v2;u0;u_correct" << endl;
	for (size_t i = 0; i < data.size(); i++)
		tecplot_file << data[i].v1 << ';' << data[i].v2 << ';' << data[i].u0 << ';' << data[i].u_correct << endl;;

	tecplot_file << "\n";
	tecplot_file.flush();
	tecplot_file.close();
}
