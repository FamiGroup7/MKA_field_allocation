#include "Mka3D.h"
#include "Mka2D_cylindrical.h"

int LosLU(double* ggl, double* ggu, double* diag, int N, int* ig, int* jg, double* f, double* q);

int main(int argc, char* argv[])
{
	setlocale(LC_ALL, "rus");

	double sourcePower = 1;
	Mka2D_cylindrical task2D("resources_field_allocation/2D/", false, true, false, true);
	task2D.power = sourcePower;
	task2D.startFull2dProcess();

	//Mka3D task("resources3D/", true, false, true, true,
	//	true, true, true);
	//task.startFullProcess();

	Mka3D task("resources_field_allocation/3D/", false, false, true, true,
		true, true, true);
	task.buildNet("sreda.txt", "sourceLocate.txt");
	task.build_xyz_nvtr_portratin_Abqx();
	double * q0 = new double[task.countRegularNodes];
	ofstream log("resources_field_allocation/log.txt");
	//task.generateGlobalMatrix();
	for (int i = 0; i < task.countRegularNodes; i++)
	{
		Point target = task.xyz_points[i];
		Point_cylindrical target2D(sqrt(target.x*target.x + target.y*target.y), target.z);
		int iKe2D = task2D.findKE(target2D);
		double solution = 0;
		if (iKe2D >= 0) {
			solution = task2D.solutionInPoint(iKe2D, target2D);
		}
		
		q0[i] = solution;
	}

	//ahTUNG!
	double lambda0 = 0.5;
	for (int ielem = 0; ielem < task.KE.size(); ielem++)
	{
		int area = task.FindAreaNumber(task.KE[ielem].uzel);
		task.CreateLocalMatrixs(ielem, lambda0 - task.sreda[area].lambda);
		task.Addition(ielem, task.di, task.ggl);
		//task.PrintLocalMatrix();
	}
	for (int i = 0; i < task.ig[task.countRegularNodes]; i++)
	{
		task.ggu[i] = task.ggl[i];
	}
	double* b_plus = new double[task.countRegularNodes];
	task.MultMatrixOnVector(q0, b_plus, task.di, task.ggl, task.ggu);

	for (int i = 0; i < task.ig[task.countRegularNodes]; i++)
	{
		task.ggu[i] = task.ggl[i] = 0;
	}
	for (int i = 0; i < task.countRegularNodes; i++) {
		task.di[i] = task.b[i] = task.q[i] = 0;
	}
	//task.upKU = task.downKU = task.leftKU = task.rightKU = task.foreKU = task.behindKU = 4;
	task.generateGlobalMatrix(lambda0);
	LosLU(task.ggl, task.ggu, task.di, task.countRegularNodes, task.ig, task.jg, b_plus, task.q);

	//solve 3D
	Mka3D taskCheck("resources_field_allocation/3D_check/", false, false, true, true,
		true, true, true);
	taskCheck.power = sourcePower;
	taskCheck.startFullProcess();
	//taskCheck.buildNet("sreda.txt", "sourceLocate.txt");
	//taskCheck.power = 1;
	////taskCheck.netOptimization();

	//taskCheck.build_xyz_nvtr_portratin_Abqx();
	//taskCheck.generateGlobalMatrix();
	//LosLU(taskCheck.ggl, taskCheck.ggu, taskCheck.di, taskCheck.countRegularNodes, 
	//	taskCheck.ig, taskCheck.jg, taskCheck.b, taskCheck.q);
	//ofstream outputCheck("resources_field_allocation/3D_check/solution.txt");
	//taskCheck.calcPogreshnost(outputCheck);


	double * q_check = new double[task.countRegularNodes];
	for (int i = 0; i < task.countRegularNodes; i++)
	{
		Point target = task.xyz_points[i];
		int iKe_checkNet = taskCheck.findKE(target);
		double solution = 0;
		if (iKe_checkNet >= 0) {
			solution = taskCheck.solutionInPoint(iKe_checkNet, target);
		}

		q_check[i] = solution;
	}

	ofstream outFIleCsv("resources_field_allocation/res.csv");
	log << setw(20) << "q0" << setw(20) << "task.q+" << setw(20) << "q_sum" << setw(20) << "q3D" << endl;
	outFIleCsv << "q0" << ';' << "task.q+" << "q_sum" << ';' << "q3D" << endl;
	double diff = 0;
	for (int i = 0; i < task.countRegularNodes; i++) {
		double sum = q0[i] + task.q[i];
		log << setw(20) << q0[i] << setw(20) << task.q[i] << setw(20) << sum << setw(20) << q_check[i] << endl;
		outFIleCsv << q0[i] << ';' << task.q[i] << ';' << sum << ';' << q_check[i] << endl;
		diff += (sum + q_check[i])*(sum + q_check[i]);
	}
	log <<"Diff="<< diff << endl;
}
