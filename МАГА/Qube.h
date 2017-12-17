#pragma once
#include "MkaUtils.h"

class Qube
{
public:
	int i_startX, i_startY, i_startZ;
	int i_nextX, i_nextY, i_nextZ;
	Qube();
	Qube(int startX, int startY, int startZ);
	Qube(int startX, int startY, int startZ, int nextX, int nextY, int nextZ);
	~Qube();
	double getWidth(double*xNet);
	double getDepth(double*yNet);
	double getHeight(double*zNet);

};
