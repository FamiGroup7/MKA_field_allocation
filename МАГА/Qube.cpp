#include "Qube.h"

Qube::Qube()
{
}

Qube::Qube(int startX, int startY, int startZ)
{
	this->i_startX = startX;
	this->i_startY = startY;
	this->i_startZ = startZ;

	this->i_nextX = -1;
	this->i_nextY = -1;
	this->i_nextZ = -1;
}

Qube::Qube(int startX, int startY, int startZ, int nextX, int nextY, int nextZ)
{
	this->i_startX = startX;
	this->i_startY = startY;
	this->i_startZ = startZ;

	this->i_nextX = nextX;
	this->i_nextY = nextY;
	this->i_nextZ = nextZ;
}

Qube::~Qube()
{
}


double Qube::getWidth(double*xNet) {
	return fabs(xNet[i_nextX] - xNet[i_startX]);
}
double Qube::getDepth(double*yNet) {
	return fabs(yNet[i_nextY] - yNet[i_startY]);
}
double Qube::getHeight(double*zNet) {
	return fabs(zNet[i_nextZ] - zNet[i_startZ]);
}

Point Qube::getCenter(double*xNet, double*yNet, double*zNet)
{
	double xCenter, yCenter, zCenter;
	xCenter = xNet[i_startX] + xNet[i_nextX];
	yCenter = yNet[i_startY] + yNet[i_nextY];
	zCenter = zNet[i_startZ] + zNet[i_nextZ];
	return Point(xCenter / 2, yCenter / 2, zCenter / 2);
}
