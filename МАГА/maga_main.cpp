#include "Mka3D.h"
#include "Mka2D_cylindrical.h"

int main(int argc, char* argv[])
{
	setlocale(LC_ALL, "rus");

	//Mka2D_cylindrical task2D("resources_field_allocation/2D/", true, true, true, true);

	Mka3D task("resources3D/", true, false, true, true,
		true, true, true);
	task.startFullProcess();

	//task.buildNet("sreda.txt", "sourceLocate.txt");
	//task.build_xyz_nvtr_portratin_Abqx();
	//double * q0 = new double[task.countRegularNodes];
	//ofstream log("resources_field_allocation/log.txt");
	//task.generateGlobalMatrix();
	//for (int i = 0; i < task.countRegularNodes; i++)
	//{
	//	Point target = task.xyz_points[i];
	//	Point_cylindrical target2D(sqrt(target.x*target.x + target.y*target.y), target.z);
	//	int iKe2D = task2D.findKE(target2D);
	//	double solution = task2D.solutionInPoint(iKe2D, target2D);
	//	log << solution << endl;
	//	q0[i] = solution;
	//}

}
