#include "Mka3D.h"
#include "Mka2D_cylindrical.h"

int main(int argc, char* argv[])
{
	setlocale(LC_ALL, "rus");
	Mka3D task(true, false, true, true,
		true, true, true);
	task.solve();

	//Mka2D_cylindrical task;
	//Mka3D task(true, "sreda.txt");
}
