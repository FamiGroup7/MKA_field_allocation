#include "Mka3D.h"
#include "Mka2D_cylindrical.h"

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

void main() {
		Mka2D_cylindrical task2D("resources2D/", true, true, true, true);
		task2D.power = 0;
		task2D.startFull2dProcess();
}

//int main(int argc, char* argv[])
//{
//	auto start = std::chrono::system_clock::now();
//	setlocale(LC_ALL, "rus");
//	ofstream log("resources_field_allocation/log.txt");
//
//	//ahTUNG!
//	double lambda0 = 1;
//	double sourcePower = 0.01;
//	Mka2D_cylindrical task2D("resources_field_allocation/2D/", true, false, false, true);
//	Mka3D task("resources_field_allocation/3D/", false, false, true, false,
//		true, true, true);
//	Mka3D taskCheck("resources_field_allocation/3D_check/", false, false, true, false,
//		true, true, true);
//	task2D.power = sourcePower;
//	//task.power = sourcePower;
//	taskCheck.power = sourcePower;
//	task2D.startFull2dProcess();
//
//	task.buildNet("sreda.txt", "sourceLocate.txt");
//	task.build_xyz_nvtr_portratin_Abqx();
//
//	double * q0 = map2dSolutionToNet(task2D, task);
//
//	//create G_delta_lambda
//	for (int ielem = 0; ielem < task.KE.size(); ielem++)
//	{
//		int area = task.FindAreaNumber(task.KE[ielem].uzel);
//		double val = 0;
//		//val = lambda0 - task.sreda[area].lambda;
//		if (area == 0) {
//			val = 0.99;
//		}
//		if (MkaUtils::equals(val, 0)) {
//			continue;
//		}
//		task.CreateLocalMatrixs(ielem, val);
//		task.Addition(ielem, task.di, task.ggl);
//	}
//	//log << "ggl" << endl;
//	for (int i = 0; i < task.ig[task.countRegularNodes]; i++)
//	{
//		//log << task.ggl[i] << endl;
//		task.ggu[i] = task.ggl[i];
//	}
//	//log << endl << "--------" << endl;
//
//	double* b_plus = new double[task.countRegularNodes];
//	memset(b_plus, 0, sizeof(double)*task.countRegularNodes);
//	task.MultMatrixOnVector(q0, b_plus, task.di, task.ggl, task.ggu);
//
//	//log << "b_plus		di" << endl;
//	//for (int i = 0; i < task.countRegularNodes; i++)
//	//{
//	//	log << b_plus[i] << '\t' << task.di[i] << endl;
//	//}
//	//log << endl << "--------" << endl;
//
//	for (int i = 0; i < task.ig[task.countRegularNodes]; i++)
//	{
//		task.ggu[i] = task.ggl[i] = 0;
//	}
//	for (int i = 0; i < task.countRegularNodes; i++) {
//		task.di[i] = task.b[i] = task.q[i] = 0;
//	}
//	task.generateGlobalMatrix(-1);
//	LosLU(task.ggl, task.ggu, task.di, task.countRegularNodes, task.ig, task.jg, b_plus, task.q);
//
//	//solve 3D
//	taskCheck.startFullProcess();
//
//	double * q_check = new double[task.countRegularNodes];
//	for (int i = 0; i < task.countRegularNodes; i++)
//	{
//		Point target = task.xyz_points[i];
//		int iKe_checkNet = taskCheck.findKE(target);
//		double solution = 0;
//		if (iKe_checkNet >= 0) {
//			solution = taskCheck.solutionInPoint(iKe_checkNet, target);
//		}
//
//		q_check[i] = solution;
//	}
//
//	ofstream outFIleCsv("resources_field_allocation/res.csv");
//	ofstream outSolutionOnUpFace("resources_field_allocation/solutionOnUpFace.csv");
//	log << setw(20) << "q0" << setw(20) << "task.q+" << setw(20) << "q_sum" << setw(20) << "q3D" << setw(20) << "diff" << endl;
//	outFIleCsv << "q0" << ';' << "task.q+" << "q_sum" << ';' << "q3D" << ';' << "diff" << endl;
//	outSolutionOnUpFace << "x" << ';' << "q0" << ';' << "task.q+" << ';' << "q_sum" << ';' << "q3D" << endl;
//	double diff = 0;
//	char buff[100];
//	for (int i = 0; i < task.countRegularNodes; i++) {
//		double sum = q0[i] + task.q[i];
//		log << setw(20) << q0[i] << setw(20) << task.q[i] << setw(20) << sum << setw(20) << q_check[i] << setw(20) << q_check[i] - sum << endl;
//		outFIleCsv << q0[i] << ';' << task.q[i] << ';' << sum << ';' << q_check[i] << ';' << q_check[i] - sum << endl;
//		diff += (sum - q_check[i])*(sum - q_check[i]);
//		if (MkaUtils::equals(task.xyz_points[i].y, 0) && MkaUtils::equals(task.xyz_points[i].z, 0)) {
//			snprintf(buff, sizeof(buff), "%15.7e;%15.7e;%15.7e;%15.7e;%15.7e", task.xyz_points[i].x, q0[i], task.q[i], sum, q_check[i]);
//			outSolutionOnUpFace << buff << endl;
//		}
//	}
//	log <<"Diff="<< diff << endl;
//	log << setw(40) << std::left << "execution time" <<
//		MkaUtils::formattingTime(std::chrono::system_clock::time_point(std::chrono::system_clock::now() - start)) << endl;
//}
