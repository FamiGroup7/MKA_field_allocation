#include "Mka3D.h"
#include "Mka2D_cylindrical.h"

ofstream logger("resources_field_allocation/log.txt");

int LosLU(double* ggl, double* ggu, double* diag, int N, int* ig, int* jg, double* f, double* q);

struct data_elem
{
	double v1, v2;
	double u_correct, u0, u_plus;
};

void output_slice(const string & tecplot_filename, Mka3D &taskCheck, double*q0, Mka3D &fieldAllocation,
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
			int iKE_allocation = fieldAllocation.findKE(p);
			if (iKE_allocation < 0) {
				throw new exception();
			}
			tmp_elem.u_plus = fieldAllocation.solutionInPoint(iKE_allocation, p);
			tmp_elem.u0 = fieldAllocation.solutionInPoint(iKE_allocation, p, q0);

			int iKE = taskCheck.findKE(p);
			if (iKE < 0) {
				throw new exception();
			}
			tmp_elem.u_correct = taskCheck.solutionInPoint(iKE, p);

			data[i * num_var_2 + j] = tmp_elem;
		}
	}

	tecplot_file << "v1;v2;u0;u_plus;u_correct" << endl;
	for (size_t i = 0; i < data.size(); i++)
		tecplot_file << data[i].v1 << ';' << data[i].v2 << ';' << data[i].u0 << ';' << data[i].u_plus << ';' << data[i].u_correct << endl;;

	tecplot_file << "\n";
	tecplot_file.flush();
	tecplot_file.close();
}

//--------------------------------------

void output_slice(const string & tecplot_filename, Mka3D &task, double*q0,
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
			int iKE = task.findKE(p);
			if (iKE < 0) {
				throw new exception();
			}
			tmp_elem.u_correct = task.solutionInPoint(iKE, p);
			tmp_elem.u0 = task.solutionInPoint(iKE, p, q0);
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

double*map2dSolutionToNet(Mka2D_cylindrical&from, Mka3D&to) {
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

//void main() {
//	Mka2D_cylindrical task2D("resources_field_allocation/2D/", false, true, false, true);
//	task2D.startFull2dProcess();
//}

//void main() {
//	Mka2D_cylindrical task2D("resources_field_allocation/2D/", false, true, false, true);
//	task2D.startFull2dProcess();
//	Mka3D task("resources_field_allocation/3D/", false, true, true, false,
//		true, true, true);
//	task.startFullProcess();
//
//	//output_slice("output/slice.csv", task, task2D, 'Z', -700, 'X', 1000, 1700.0, 200, 'Y', 500.0, 500.0, 1);
//		output_slice("output/sliceZ-450_Y0.csv", task, task2D,  'Z', -450, 'X', 0, 600.0, 150, 'Y', 0.0, 0.0, 1);
//		output_slice("output/sliceZ-450_Y450.csv", task, task2D,  'Z', -450, 'X', 0, 600.0, 150, 'Y', 450.0, 450.0, 1);
//}

//void main() {
//	ofstream log("resources_field_allocation/log.txt");
//	double sourcePower = 1;
//	Mka2D_cylindrical task2D("resources_field_allocation/2D/", false, true, false, true);
//	Mka3D taskCheck("resources_field_allocation/3D/", false, false, true, false,
//		true, true, true);
//	task2D.startFull2dProcess();
//	//solve 3D
//	taskCheck.startFullProcess();
//	//output_slice("outputs/slice3d.txt", taskCheck, 'Y', -220.0, 'X', -500.0, 500.0, 200, 'Z', -500.0, 500.0, 200);
//
//	double*q0 = map2dSolutionToNet(task2D, taskCheck);
//	double*q_check = taskCheck.q;
//
//	ofstream outFIleCsv("resources_field_allocation/res.csv");
//	ofstream outSolutionOnUpFace("resources_field_allocation/solutionOnUpFace.csv");
//	log << setw(20) << "q0" << setw(20) << "q3D" << setw(20) << "diff" << setw(20) << "x" << setw(20) << "y" << setw(20) << "z" << endl;
//	outFIleCsv << "q0" << ';' << "q3D" << ';' << "diff" << endl;
//	outSolutionOnUpFace << "x" << ';' << "q0" << ';' << "q3D" << endl;
//	double diff = 0, accurate = 0;
//	char buff[100];
//	for (int i = 0; i < taskCheck.countRegularNodes; i++) {
//		accurate += q_check[i] * q_check[i];
//		double sum = q0[i];
//		log << setw(20) << q0[i]  << setw(20) << q_check[i] << setw(20) << q_check[i] - sum << setw(20) << taskCheck.xyz_points[i].x << setw(20) << taskCheck.xyz_points[i].y << setw(20) << taskCheck.xyz_points[i].z << endl;
//		outFIleCsv << q0[i] << ';' << q_check[i] << ';' << q_check[i] - sum << endl;
//		diff += (sum - q_check[i])*(sum - q_check[i]);
//		if (MkaUtils::equals(taskCheck.xyz_points[i].y, 0) && MkaUtils::equals(taskCheck.xyz_points[i].z, 0)) {
//			snprintf(buff, sizeof(buff), "%e;%e;%e", taskCheck.xyz_points[i].x, q0[i], q_check[i]);
//			outSolutionOnUpFace << buff << endl;
//		}
//	}
//	log << "Diff=" << diff << endl;
//	log << "||Diff/accurate||=" << diff / accurate << endl;
//}

void generageG_lambda(double lambda0, Mka3D&simple) {
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

	logger << "ggl" << endl;
	for (int i = 0; i < simple.ig[simple.countRegularNodes]; i++)
	{
		//logger << simple.ggl[i] << endl;
	}
	logger << endl << "--------" << endl;
}

//return solution u0
void fieldAllocation(double lambda0, double*q0, Mka3D&simple) {
	generageG_lambda(lambda0, simple);

	double* b_plus = new double[simple.countRegularNodes];
	memset(b_plus, 0, sizeof(double)*simple.countRegularNodes);
	simple.MultMatrixOnVector(q0, b_plus, simple.di, simple.ggl, simple.ggu);

	logger << "b_plus		di" << endl;
	for (int i = 0; i < simple.countRegularNodes; i++)
	{
		//logger << b_plus[i] << '\t' << simple.di[i] << endl;
	}
	logger << endl << "--------" << endl;

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

////without 2d
//int main(int argc, char* argv[])
//{
//	auto start = std::chrono::system_clock::now();
//	setlocale(LC_ALL, "rus");
//
//	////ahTUNG!
//	double lambda0 = 1;
//	Mka3D task_instead2D("resources_field_allocation/3D/", false, true, true, false,
//		true, true, true);
//	task_instead2D.buildNet("sreda.txt", "sourceLocate.txt");
//	task_instead2D.build_xyz_nvtr_portratin_Abqx();
//	task_instead2D.generateGlobalMatrix(lambda0);
//	LosLU(task_instead2D.ggl, task_instead2D.ggu, task_instead2D.di, task_instead2D.countRegularNodes, 
//		task_instead2D.ig, task_instead2D.jg, task_instead2D.b, task_instead2D.q);
//	Mka3D task("resources_field_allocation/3D/", false, false, true, false,
//		true, true, true);
//	Mka3D taskCheck("resources_field_allocation/3D_check/", false, false, true, false,
//		true, true, true);
//
//	task.buildNet("sreda.txt", "sourceLocate.txt");
//	task.build_xyz_nvtr_portratin_Abqx();
//
//	double*q0 = task_instead2D.q;
//	fieldAllocation(lambda0, q0, task);
//
//	////check
//	taskCheck.startFullProcess();
//	double * q_check = new double[task.countRegularNodes];
//	for (int i = 0; i < task.countRegularNodes; i++)
//	{
//		Point target = task.xyz_points[i];
//		int iKe_checkNet = taskCheck.findKE(target);
//		double solution = 0;
//		if (iKe_checkNet >= 0) {
//			solution = taskCheck.solutionInPoint(iKe_checkNet, target);
//		}
//		else {
//			i = i;
//			int iKe = task.findKE(target);
//			if (iKe < 0) {
//				i = i;
//			}
//			else {
//				double val = task.solutionInPoint(iKe, target);
//				i = i;
//			}
//		}
//
//		q_check[i] = solution;
//	}
//
//	ofstream outFIleCsv("resources_field_allocation/res.csv");
//	ofstream outSolutionOnUpFace("resources_field_allocation/solutionOnUpFace.csv");
//	logger << setw(20) << "q0" << setw(20) << "task.q+" << setw(20) << "q_sum" << setw(20) << "q3D" << setw(20) << "diff" << setw(20) << "x" << setw(20) << "y" << setw(20) << "z" << endl;
//	outFIleCsv << "q0" << ';' << "task.q+" << "q_sum" << ';' << "q3D" << ';' << "diff" << endl;
//	outSolutionOnUpFace << "x" << ';' << "q0" << ';' << "task.q+" << ';' << "q_sum" << ';' << "q3D" << endl;
//	double diff = 0, accurate = 0;
//	char buff[100];
//	for (int i = 0; i < task.countRegularNodes; i++) {
//		accurate += q_check[i] * q_check[i];
//		double sum = q0[i] + task.q[i];
//		logger << setw(20) << q0[i] << setw(20) << task.q[i] << setw(20) << sum << setw(20) << q_check[i] << setw(20) << q_check[i] - sum << setw(20) << task.xyz_points[i].x << setw(20) << task.xyz_points[i].y << setw(20) << task.xyz_points[i].z << endl;
//		outFIleCsv << q0[i] << ';' << task.q[i] << ';' << sum << ';' << q_check[i] << ';' << q_check[i] - sum << endl;
//		diff += (sum - q_check[i])*(sum - q_check[i]);
//		if (MkaUtils::equals(task.xyz_points[i].y, 0) && MkaUtils::equals(task.xyz_points[i].z, 0)) {
//			snprintf(buff, sizeof(buff), "%e;%e;%e;%e;%e", task.xyz_points[i].x, q0[i], task.q[i], sum, q_check[i]);
//			outSolutionOnUpFace << buff << endl;
//		}
//	}
//	logger << "Diff=" << diff << endl;
//	logger << "||Diff/accurate||=" << diff / accurate << endl;
//	logger << setw(40) << std::left << "execution time" <<
//		MkaUtils::formattingTime(std::chrono::system_clock::time_point(std::chrono::system_clock::now() - start)) << endl;
//
//	output_slice("output/sliceZ-450_Y0.csv", taskCheck, q0, task, 'Z', -450, 'X', 0, 600.0, 150, 'Y', 0.0, 0.0, 1);
//	output_slice("output/sliceZ-450_Y450.csv", taskCheck, q0, task, 'Z', -450, 'X', 0, 600.0, 150, 'Y', 450.0, 450.0, 1);
//}

//check themselfs
int main(int argc, char* argv[])
{
	auto start = std::chrono::system_clock::now();
	setlocale(LC_ALL, "rus");

	//ahTUNG!
	double lambda0 = 1;
	Mka3D task("resources_field_allocation/3D/", false, false, true, false,
		true, true, true);
	Mka3D taskCheck("resources_field_allocation/3D/", false, false, true, false,
		true, true, true);

	task.startFullProcess();
	taskCheck.startFullProcess();

	double * q_check = new double[task.countRegularNodes];
	for (int i = 0; i < task.countRegularNodes; i++)
	{
		Point target = task.xyz_points[i];
		int iKe_checkNet = taskCheck.findKE(target);
		double solution = 0;
		if (iKe_checkNet >= 0) {
			solution = taskCheck.solutionInPoint(iKe_checkNet, target);
		}
		else {
			i = i;
			int iKe = task.findKE(target);
			if (iKe < 0) {
				i = i;
			}
			else {
				double val = task.solutionInPoint(iKe, target);
				i = i;
			}
		}

		q_check[i] = solution;
	}

	logger << setw(40) << std::left << "execution time" <<
		MkaUtils::formattingTime(std::chrono::system_clock::time_point(std::chrono::system_clock::now() - start)) << endl;

	output_slice("output/sliceZ-450_Y0.csv", taskCheck, task.q, 'Z', -450, 'X', 0, 600.0, 150, 'Y', 0.0, 0.0, 1);
	output_slice("output/sliceZ-450_Y450.csv", taskCheck, task.q, 'Z', -450, 'X', 0, 600.0, 150, 'Y', 450.0, 450.0, 1);
}

//int main(int argc, char* argv[])
//{
//	auto start = std::chrono::system_clock::now();
//	setlocale(LC_ALL, "rus");
//
//	//ahTUNG!
//	double lambda0 = 1;
//	Mka2D_cylindrical task2D("resources_field_allocation/2D/", false, true, false, true);
//	Mka3D task("resources_field_allocation/3D/", false, false, true, false,
//		true, true, true);
//	Mka3D taskCheck("resources_field_allocation/3D_check/", false, false, true, false,
//		true, true, true);
//	task2D.startFull2dProcess();
//
//	task.buildNet("sreda.txt", "sourceLocate.txt");
//	task.build_xyz_nvtr_portratin_Abqx();
//
//	double *q0 = map2dSolutionToNet(task2D, task);
//	fieldAllocation(lambda0, q0, task);
//
//	//check
//	taskCheck.startFullProcess();
//	double * q_check = new double[task.countRegularNodes];
//	for (int i = 0; i < task.countRegularNodes; i++)
//	{
//		Point target = task.xyz_points[i];
//		int iKe_checkNet = taskCheck.findKE(target);
//		double solution = 0;
//		if (iKe_checkNet >= 0) {
//			solution = taskCheck.solutionInPoint(iKe_checkNet, target);
//		}
//		else {
//			i = i;
//		}
//
//		q_check[i] = solution;
//	}
//
//	ofstream outFIleCsv("resources_field_allocation/res.csv");
//	ofstream outSolutionOnUpFace("resources_field_allocation/solutionOnUpFace.csv");
//	logger << setw(20) << "q0" << setw(20) << "task.q+" << setw(20) << "q_sum" << setw(20) << "q3D" << setw(20) << "diff" << setw(20) << "x" << setw(20) << "y" << setw(20) << "z" << endl;
//	outFIleCsv << "q0" << ';' << "task.q+" << "q_sum" << ';' << "q3D" << ';' << "diff" << endl;
//	outSolutionOnUpFace << "x" << ';' << "q0" << ';' << "task.q+" << ';' << "q_sum" << ';' << "q3D" << endl;
//	double diff = 0, accurate = 0;
//	char buff[100];
//	for (int i = 0; i < task.countRegularNodes; i++) {
//		accurate += q_check[i] * q_check[i];
//		double sum = q0[i] + task.q[i];
//		logger << setw(20) << q0[i] << setw(20) << task.q[i] << setw(20) << sum << setw(20) << q_check[i] << setw(20) << q_check[i] - sum <<setw(20)<<task.xyz_points[i].x << setw(20) << task.xyz_points[i].y << setw(20) << task.xyz_points[i].z << endl;
//		outFIleCsv << q0[i] << ';' << task.q[i] << ';' << sum << ';' << q_check[i] << ';' << q_check[i] - sum << endl;
//		diff += (sum - q_check[i])*(sum - q_check[i]);
//		if (MkaUtils::equals(task.xyz_points[i].y, 0) && MkaUtils::equals(task.xyz_points[i].z, 0)) {
//			snprintf(buff, sizeof(buff), "%e;%e;%e;%e;%e", task.xyz_points[i].x, q0[i], task.q[i], sum, q_check[i]);
//			outSolutionOnUpFace << buff << endl;
//		}
//	}
//	logger <<"Diff="<< diff << endl;
//	logger << "||Diff/accurate||=" << diff / accurate << endl;
//	logger << setw(40) << std::left << "execution time" <<
//		MkaUtils::formattingTime(std::chrono::system_clock::time_point(std::chrono::system_clock::now() - start)) << endl;
//	
//	output_slice("output/sliceZ-450_Y0.csv", taskCheck, q0, task, 'Z', -450, 'X', 0, 600.0, 150, 'Y', 0.0, 0.0, 1);
//	output_slice("output/sliceZ-450_Y450.csv", taskCheck, q0, task, 'Z', -450, 'X', 0, 600.0, 150, 'Y', 450.0, 450.0, 1);
//}
