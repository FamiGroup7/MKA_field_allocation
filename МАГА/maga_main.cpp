#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <locale>
#include <bitset>
#include <Windows.h>
#include "Edge.h"

const size_t BIT_SIZE = 13;
const size_t LEFT_UP = 0;
const size_t LEFT_DOWN = 1;
const size_t RIGHT_UP = 2;
const size_t RIGHT_DOWN = 3;
const size_t FORE_UP = 4;
const size_t FORE_DOWN = 5;
const size_t BACK_UP = 6;
const size_t BACK_DOWN = 7;
const size_t LEFT_FORE = 8;
const size_t RIGHT_FORE = 9;
const size_t LEFT_BACK = 10;
const size_t RIGHT_BACK = 11;
const size_t IS_REGULAR = 12;

const bool DEBUG = true;

using namespace std;

struct colour // цвет точки
{
	int red;
	int green;
	int blue;
};

struct point // точка
{
	double x, y, z;

	point(double x, double y, double z)
		: x(x),
		  y(y),
		  z(z)
	{
	}

	point(): x(0), y(0), z(0)
	{
	}

	friend bool operator==(const point& lhs, const point& rhs)
	{
		return MkaUtils::equals(lhs.x, rhs.x)
			&& MkaUtils::equals(lhs.y, rhs.y)
			&& MkaUtils::equals(lhs.z, rhs.z);
	}

	friend bool operator!=(const point& lhs, const point& rhs)
	{
		return !(lhs == rhs);
	}
};

struct locateOfPoint
{
	int i, j, k;

	locateOfPoint(int i, int j, int k)
		: i(i),
		  j(j),
		  k(k)
	{
	}

	locateOfPoint(): i(0), j(0), k(0)
	{
	}
};

struct field
{
	double x1, x2, y1, y2, z1, z2, lambda, gamma;

	field(double x1, double x2, double y1, double y2, double z1, double z2, double lambda, double gamma)
		: x1(x1),
		  x2(x2),
		  y1(y1),
		  y2(y2),
		  z1(z1),
		  z2(z2),
		  lambda(lambda),
		  gamma(gamma)
	{
	}

	field()
	{
	}
};

struct nvtr
{
	int uzel[8], numberField;

	nvtr(int number_field, int newNodes[8])
		: numberField(number_field)
	{
		memcpy(uzel, newNodes, 8);
	}

	nvtr()
	{
	}
};

struct neighbor {
	int index;
	double weight;

	neighbor(int ind, double w) :index(ind), weight(w) {
	}

	friend bool operator< (const neighbor &left, const neighbor &right) {
		return left.index < right.index;
	}
};

struct sigmStruct3D {
	int terminalNode;
	set<neighbor> neighbors;
} *tmpSigm;
vector<sigmStruct3D> sigmNewT;

int *igT, *jgT;
double* ggT;

int LosLU(double* ggl, double* ggu, double* diag, int N, int* ig, int* jg, double* f, double* q);

vector<point> xyz_points;
vector<nvtr> KE;
vector<field> sreda;
map<Point, Edge> termNodeOnEdge;


double leftX, rightX, leftY, rightY, leftZ, rightZ;
double koordSourceX, koordSourceY, koordSourceZ;
double *xNet, *yNet, *zNet;
int nX, nY, nZ;
char*** matrixNode;
bitset<BIT_SIZE>*** newNodes;
int nColT, kolvoRegularNode;

int *ig, *jg;
double *di, *b, *q, *ggl, *ggu;
double localB[8], localMatrix[8][8];
double G1[2][2] = {{1,-1},{-1,1}};
double M1[2][2] = {{1 / 3.,1 / 6.},{1 / 6.,1 / 3.}};
double M2[4][4] = {
	{4,2,2,1},
	{2,4,1,2},
	{2,1,4,2},
	{1,2,2,4}
};
ofstream output("solution.txt");
ofstream logger("log.txt");

void inputConfig()
{
}

void logError(char* message) {
	cout << message << endl;
	throw exception(message);
}

void GenerateNetLikeTelma(set<double>& mas, ifstream& fileNet)
{
	double firstElement, startPosition, endPosition, position, lastPosition, stepReal;;
	int i, numberOfInterval;
	fileNet >> firstElement >> numberOfInterval;
	double* intervals = new double[numberOfInterval + 1];
	double* sizeOfSteps = new double[numberOfInterval];
	double* mnojiteli = new double[numberOfInterval];
	int* napravlenie = new int[numberOfInterval];
	intervals[0] = firstElement;
	for (i = 1; i < numberOfInterval + 1; i++)
	{
		fileNet >> intervals[i];
	}
	for (i = 0; i < numberOfInterval; i++)
	{
		fileNet >> sizeOfSteps[i];
	}
	for (i = 0; i < numberOfInterval; i++)
	{
		fileNet >> mnojiteli[i];
		if (mnojiteli[i] < 1)
		{
			cout << "Ошибка описания сетки по оси. Коэфициент разрядки должен быть больше или равен 1";
			system("pause");
			exit(1);
		}
	}
	for (i = 0; i < numberOfInterval; i++)
	{
		fileNet >> napravlenie[i];
	}
	for (i = 0; i < numberOfInterval; i++)
	{
		if (napravlenie[i] == -1)
		{
			startPosition = intervals[i + 1];
			endPosition = intervals[i];
		}
		else
		{
			startPosition = intervals[i];
			endPosition = intervals[i + 1];
		}
		lastPosition = startPosition;
		position = startPosition + napravlenie[i] * sizeOfSteps[i];
		while (position * napravlenie[i] < endPosition * napravlenie[i])
		{
			mas.insert(position);
			stepReal = fabs(lastPosition - position) * mnojiteli[i];
			lastPosition = position;
			position += napravlenie[i] * stepReal;
		}
		/*if (fabs(lastPosition - endPosition) < sizeOfSteps[i]*0.9 && lastPosition != startPosition)
		{
			mas.erase(mas.find(lastPosition));
		}*/
	}
}

int indexXYZ(point goal)
{
	for (int i = 0; i < xyz_points.size(); i++)
	{
		if (goal == xyz_points[i])
			return i;
	}

	cout << "Ошибка. Не найдена точка (" << goal.x << "," << goal.y << "," << goal.z << ")" << endl;
	throw new exception();
	system("pause");
	exit(1);
}

int indexXYZ(double x, double y, double z)
{
	point goal(x, y, z);
	return indexXYZ(goal);
}

double LikeASquare(double x, double y)
{
	if (x > y)
		return (1 - y / x);
	else
		return (1 - x / y);
}

double otn(double val1, double val2) {
	if (val1 > val2)return val2 / val1;
	return val1 / val2;
}

double LikeACube(double x, double y, double z) {
	return (otn(x, y) + otn(x, z) + otn(y, z)) / 3.0;
}

int FindLocate(double* massiv, int razm, double x)
{
	int i;
	for (i = 0; i < razm; i++)
	{
		if (massiv[i] == x)
			return i;
	}
	return -1;
}

bool hasLeft(bitset<BIT_SIZE> node) {
	return node.test(LEFT_BACK) || node.test(LEFT_FORE)
		|| node.test(LEFT_DOWN) || node.test(LEFT_UP);
}

bool hasRight(bitset<BIT_SIZE> node) {
	return node.test(RIGHT_BACK) || node.test(RIGHT_FORE)
		|| node.test(RIGHT_DOWN) || node.test(RIGHT_UP);
}

bool hasUp(bitset<BIT_SIZE> node) {
	return node.test(LEFT_UP) || node.test(RIGHT_UP)
		|| node.test(BACK_UP) || node.test(FORE_UP);
}

bool hasDown(bitset<BIT_SIZE> node) {
	return node.test(LEFT_DOWN) || node.test(RIGHT_DOWN)
		|| node.test(BACK_DOWN) || node.test(FORE_DOWN);
}

bool hasBack(bitset<BIT_SIZE> node) {
	return node.test(LEFT_BACK) || node.test(RIGHT_BACK)
		|| node.test(BACK_UP) || node.test(BACK_DOWN);
}

bool hasFore(bitset<BIT_SIZE> node) {
	return node.test(LEFT_FORE) || node.test(RIGHT_FORE)
		|| node.test(FORE_UP) || node.test(FORE_DOWN);
}

bool hasXLines(bitset<BIT_SIZE> node) {
	return hasLeft(node) && hasRight(node);
}

bool hasYLines(bitset<BIT_SIZE> node) {
	return hasBack(node) && hasFore(node);
}

bool hasZLines(bitset<BIT_SIZE> node) {
	return hasUp(node) && hasDown(node);
}

bool hasAll(bitset<BIT_SIZE> node) {
	return hasXLines(node) && hasYLines(node) && hasZLines(node);
}

bool canOptimize(int directionX, int directionY, int directionZ, bitset<BIT_SIZE> node) {
	if (directionX == 1)
		if (directionY == 1)
			if (directionZ == 1)
				return node.test(RIGHT_BACK) && node.test(RIGHT_UP) && node.test(BACK_UP);
			else
				return node.test(RIGHT_BACK) && node.test(RIGHT_DOWN) && node.test(BACK_DOWN);
		else
			if (directionZ == 1)
				return node.test(RIGHT_FORE) && node.test(RIGHT_UP) && node.test(FORE_UP);
			else
				return node.test(RIGHT_FORE) && node.test(RIGHT_DOWN) && node.test(FORE_DOWN);
	else
		if (directionY == 1)
			if (directionZ == 1)
				return node.test(LEFT_BACK) && node.test(LEFT_UP) && node.test(BACK_UP);
			else
				return node.test(LEFT_BACK) && node.test(LEFT_DOWN) && node.test(BACK_DOWN);
		else
			if (directionZ == 1)
				return node.test(LEFT_FORE) && node.test(LEFT_UP) && node.test(FORE_UP);
			else
				return node.test(LEFT_FORE) && node.test(LEFT_DOWN) && node.test(FORE_DOWN);
	logError("something goes wrong in \"canOptimize\"");
}

int findNextX(int directionX, int directionY, int directionZ,
	int startPositionX, int endPositionX, int posJ, int posT) {
	int posK;
	for (int k = startPositionX + 1; k <= endPositionX * directionX; k++)
	{
		posK = k * directionX;
		if (canOptimize(-1 * directionX, directionY, directionZ, newNodes[posK][posJ][posT]))
		{
			return posK;
		}
	}
	return -1;
}

int findNextY(int directionX, int directionY, int directionZ,
	int startPositionY, int endPositionY, int posI, int posT) {
	int posK;
	for (int k = startPositionY + 1; k <= endPositionY * directionY; k++)
	{
		posK = k * directionY;
		if (canOptimize(directionX, -1 * directionY, directionZ, newNodes[posI][posK][posT]))
		{
			return posK;
		}
	}
	return -1;
}

int findNextZ(int directionX, int directionY, int directionZ,
	int startPositionZ, int endPositionZ, int posI, int posJ) {
	int posK;
	for (int k = startPositionZ + 1; k <= endPositionZ * directionZ; k++)
	{
		posK = k * directionX;
		if (canOptimize(directionX, directionY, -1 * directionZ, newNodes[posI][posJ][posK]))
		{
			return posK;
		}
	}
	return -1;
}

void deletePlaneX(int xPlane, int y1, int y2, int z1, int z2) {
	logger << "Deleted plane x = " << xPlane << ", y1=" << y1 << ", y2=" << y2 << ", z1=" << z1 << ", z2=" << z2 << endl;
	
	newNodes[xPlane][y1][z1].reset(IS_REGULAR);
	newNodes[xPlane][y1][z2].reset(IS_REGULAR);
	newNodes[xPlane][y2][z1].reset(IS_REGULAR);
	newNodes[xPlane][y2][z2].reset(IS_REGULAR);

	newNodes[xPlane][y1][z1].reset(BACK_UP);
	newNodes[xPlane][y1][z2].reset(BACK_DOWN);
	newNodes[xPlane][y2][z1].reset(FORE_UP);
	newNodes[xPlane][y2][z2].reset(FORE_DOWN);

	//нижняя
	if (!newNodes[xPlane][y1][z1].test(BACK_DOWN))
		newNodes[xPlane][y1][z1].reset(LEFT_BACK).reset(RIGHT_BACK);
	if (!newNodes[xPlane][y2][z1].test(FORE_DOWN))
		newNodes[xPlane][y2][z1].reset(LEFT_FORE).reset(RIGHT_FORE);
	//верхняя
	if (!newNodes[xPlane][y1][z2].test(BACK_UP))
		newNodes[xPlane][y1][z2].reset(LEFT_BACK).reset(RIGHT_BACK);
	if (!newNodes[xPlane][y2][z2].test(FORE_UP))
		newNodes[xPlane][y2][z2].reset(LEFT_FORE).reset(RIGHT_FORE);
	//ближняя
	if (!newNodes[xPlane][y1][z1].test(FORE_UP))
		newNodes[xPlane][y1][z1].reset(LEFT_UP).reset(RIGHT_UP);
	if (!newNodes[xPlane][y1][z2].test(FORE_DOWN))
		newNodes[xPlane][y1][z2].reset(LEFT_DOWN).reset(RIGHT_DOWN);
	//дальняя
	if (!newNodes[xPlane][y2][z1].test(BACK_UP))
		newNodes[xPlane][y2][z1].reset(LEFT_UP).reset(RIGHT_UP);
	if (!newNodes[xPlane][y2][z2].test(BACK_DOWN))
		newNodes[xPlane][y2][z2].reset(LEFT_DOWN).reset(RIGHT_DOWN);
}

void addTermNodeToMap(int middleX, int x1, int x2, int y, int z) {
	termNodeOnEdge.insert(pair<Point, Edge>(
		Point(xNet[middleX], yNet[y], zNet[z]),
		Edge(
			Point(xNet[x1], yNet[y], zNet[z]),
			Point(xNet[x2], yNet[y], zNet[z]))
		));
}

void OptimizationQuarterX(int directionX, int directionY, int directionZ,
	int startX, int startY, int startZ, int endX, int endY, int endZ)
{
	int i, j, t, k, granX, granY, posI, posJ, posT, posK;
	double width_1, width_2;
	bool error;
	for (t = startZ * directionZ; t < endZ * directionZ; t++)
	{
		posT = t * directionZ;
		for (j = startY * directionY; j < endY * directionY; j++)
		{
			posJ = j * directionY;
			for (i = startX * directionX; i < endX * directionX - 1; i++)
			{
				posI = i * directionX;
				//обработка вершины при попытке расширения
				if (canOptimize(directionX, directionY, directionZ, newNodes[posI][posJ][posT]))
				{
					width_1 = width_2 = granX = granY =  0;

					//ширина 1
					int nextX = findNextX(directionX, directionY, directionZ, i, endX, posJ, posT);
					if (nextX == -1) {
						logError("error optimizationX: no width in KE "/* + posI + " " + posJ + " " + posT*/);
					}
					width_1 = fabs(xNet[nextX] - xNet[posI]);

					//глубина 1
					int nextY = findNextY(directionX, directionY, directionZ, j, endY, posI, posT);
					//глубина 2
					int checkY = findNextY(directionX, directionY, directionZ, j, endY, nextX*directionX, posT);

					//высота 1
					int nextZ = findNextZ(directionX, directionY, directionZ, t, endZ, posI, posJ);
					//высота 2
					int checkZ = findNextZ(directionX, directionY, directionZ, t, endZ, nextX*directionX, posJ);

					if (nextY == checkY && nextZ == checkZ) {
						double deep = fabs(yNet[nextY] - yNet[posJ]);
						double height = fabs(zNet[nextZ] - zNet[posT]);
						//ширина 2
						int nextX2 = findNextX(directionX, directionY, directionZ, nextX*directionX, endX, posJ, posT);
						width_2 = fabs(xNet[nextX2] - xNet[posI]);

						if (LikeACube(width_1, deep, height) < LikeACube(width_2, deep, height)) {
							if (directionY == 1)
								if (directionZ == 1)
									deletePlaneX(nextX, posJ, nextY, posT, nextZ);
								else
									deletePlaneX(nextX, posJ, nextY, nextZ, posT);
							else
								if (directionZ == 1)
									deletePlaneX(nextX, nextY, posJ, posT, nextZ);
								else
									deletePlaneX(nextX, nextY, posJ, nextZ, posT);

							addTermNodeToMap(nextX, posI, nextX2, posJ, posT);
							addTermNodeToMap(nextX, posI, nextX2, posJ, nextZ);
							addTermNodeToMap(nextX, posI, nextX2, nextY, posT);
							addTermNodeToMap(nextX, posI, nextX2, nextY, nextZ);
						}
					}
				}
			}
		}
	}
}

void initNet(double* xNet, int nX, double* yNet, int nY, double* zNet, int nZ) {
	newNodes = new bitset<BIT_SIZE>**[nX];
	for (int i = 0; i < nX; i++)
	{
		newNodes[i] = new bitset<BIT_SIZE>*[nY];
		for (int j = 0; j < nY; j++)
		{
			newNodes[i][j] = new bitset<BIT_SIZE>[nZ];
			for (int k = 0; k < nZ; k++)
			{
				bitset<BIT_SIZE> tmp;
				tmp.set();
				if (k == 0) {
					tmp.reset(LEFT_DOWN);
					tmp.reset(RIGHT_DOWN);
					tmp.reset(FORE_DOWN);
					tmp.reset(BACK_DOWN);
				}
				else if (k == nZ - 1) {
					tmp.reset(LEFT_UP);
					tmp.reset(RIGHT_UP);
					tmp.reset(FORE_UP);
					tmp.reset(BACK_UP);
				}

				if (j == 0) {
					tmp.reset(LEFT_FORE);
					tmp.reset(RIGHT_FORE);
					tmp.reset(FORE_UP);
					tmp.reset(FORE_DOWN);
				}
				else if (j == nY - 1) {
					tmp.reset(LEFT_BACK);
					tmp.reset(RIGHT_BACK);
					tmp.reset(BACK_UP);
					tmp.reset(BACK_DOWN);
				}

				if (i == 0) {
					tmp.reset(LEFT_FORE);
					tmp.reset(LEFT_BACK);
					tmp.reset(LEFT_UP);
					tmp.reset(LEFT_DOWN);
				}
				else if (i == nX - 1) {
					tmp.reset(RIGHT_FORE);
					tmp.reset(RIGHT_BACK);
					tmp.reset(RIGHT_UP);
					tmp.reset(RIGHT_DOWN);
				}
				newNodes[i][j][k] = tmp;
			}
		}
	}
}

void DivideArea(double* xNet, int nX, double* yNet, int nY, double* zNet, int nZ)
{
	double height_1, height_2, width_1, width_2;
	int locateSourceX, locateSourceY, locateSourceZ;
	for (int indSreda = 0; indSreda < sreda.size(); indSreda++)
	{
		if (koordSourceX <= sreda[indSreda].x1)
			locateSourceX = FindLocate(xNet, nX, sreda[indSreda].x1);
		else if (koordSourceX >= sreda[indSreda].x2)
			locateSourceX = FindLocate(xNet, nX, sreda[indSreda].x2);
		else
			locateSourceX = FindLocate(xNet, nX, koordSourceX);

		if (koordSourceY <= sreda[indSreda].y1)
			locateSourceY = FindLocate(yNet, nY, sreda[indSreda].y1);
		else if (koordSourceY >= sreda[indSreda].y2)
			locateSourceY = FindLocate(yNet, nY, sreda[indSreda].y2);
		else
			locateSourceY = FindLocate(yNet, nY, koordSourceY);

		if (koordSourceZ <= sreda[indSreda].z1)
			locateSourceZ = FindLocate(zNet, nZ, sreda[indSreda].z1);
		else if (koordSourceZ >= sreda[indSreda].z2)
			locateSourceZ = FindLocate(zNet, nZ, sreda[indSreda].z2);
		else
			locateSourceY = FindLocate(zNet, nZ, koordSourceZ);

		OptimizationQuarterX(1, 1, 1, locateSourceX, locateSourceY, locateSourceZ, FindLocate(xNet, nX, sreda[indSreda].x2), 
			FindLocate(yNet, nY, sreda[indSreda].y2), FindLocate(zNet, nZ, sreda[indSreda].z2));
	}
	//duplicatingOfXY();
}

void inputNet()
{
	int i, numberFields;
	field tempField;
	set<double> xTemp, yTemp, zTemp;
	set<double>::const_iterator it;
	ifstream inpSreda("sreda.txt");
	inpSreda >> numberFields;
	double areaOfSreda = 0;
	for (i = 0; i < numberFields; i++)
	{
		inpSreda >> tempField.x1 >> tempField.x2 >> tempField.y1 >> tempField.y2 >> tempField.z1 >> tempField.z2 >> tempField.lambda >> tempField.gamma;
		sreda.push_back(tempField);
		xTemp.insert(tempField.x1);
		xTemp.insert(tempField.x2);
		yTemp.insert(tempField.y1);
		yTemp.insert(tempField.y2);
		zTemp.insert(tempField.z1);
		zTemp.insert(tempField.z2);
		areaOfSreda += (tempField.x2 - tempField.x1) * (tempField.y2 - tempField.y1) * (tempField.z2 - tempField.z1);
	}
	GenerateNetLikeTelma(xTemp, inpSreda);
	GenerateNetLikeTelma(yTemp, inpSreda);
	GenerateNetLikeTelma(zTemp, inpSreda);
	set<double>::reverse_iterator rit;
	it = xTemp.begin();
	leftX = *it;
	rit = xTemp.rbegin();
	rightX = *rit;

	it = yTemp.begin();
	leftY = *it;
	rit = yTemp.rbegin();
	rightY = *rit;

	it = zTemp.begin();
	leftZ = *it;
	rit = zTemp.rbegin();
	rightZ = *rit;

	ifstream sourceFile("sourceLocate.txt");
	sourceFile >> koordSourceX >> koordSourceY >> koordSourceZ;
	if (koordSourceX > leftX && koordSourceX < rightX)
		xTemp.insert(koordSourceX);
	if (koordSourceY > leftY && koordSourceY < rightY)
		yTemp.insert(koordSourceY);
	if (koordSourceZ > leftZ && koordSourceZ < rightZ)
		zTemp.insert(koordSourceZ);

	if (fabs(areaOfSreda - (rightX - leftX) * (rightY - leftY) * (rightZ - leftZ)) > 1e-15)
	{
		cout << "Ошибка. Не правильно заданы подобласти. Общая площадь не равна сумме площадей подобластей." << endl;
		system("pause");
		exit(1);
	}

	nX = xTemp.size();
	xNet = new double[nX];
	nY = yTemp.size();
	yNet = new double[nY];
	nZ = zTemp.size();
	zNet = new double[nZ];
	i = 0;
	for (it = xTemp.begin(); it != xTemp.end(); ++it , i++)
	{
		xNet[i] = *it;
	}
	i = 0;
	for (it = yTemp.begin(); it != yTemp.end(); ++it , i++)
	{
		yNet[i] = *it;
	}
	i = 0;
	for (it = zTemp.begin(); it != zTemp.end(); ++it , i++)
	{
		zNet[i] = *it;
	}
}

void generatePortrait()
{
	int kolvoRegularNode = xyz_points.size();
	int countLocalIndex = 8;
	set<size_t>* portrait = new set<size_t>[kolvoRegularNode];
	for (size_t k = 0; k < KE.size(); k++)
	{
		for (size_t i = 0; i < countLocalIndex; i++)
		{
			size_t a = KE[k].uzel[i];
			for (size_t j = 0; j < i; j++)
			{
				size_t b = KE[k].uzel[j];
				// Если оба узла не являются терминальными
				if (a < kolvoRegularNode && b < kolvoRegularNode)
				{
					if (b > a)
						portrait[b].insert(a);
					else
						portrait[a].insert(b);
				}
				/*else if (a >= kolvoRegularNode && b < kolvoRegularNode)
				{
					for (size_t mu = igT[a - kolvoRegularNode]; mu < igT[a - kolvoRegularNode + 1]; mu++)
					{
						size_t pos_a = jgT[mu];
						if (b != pos_a)
						{
							if (b > pos_a)
								portrait[b].insert(pos_a);
							else
								portrait[pos_a].insert(b);
						}
					}
				}
				else if (a < kolvoRegularNode && b >= kolvoRegularNode)
				{
					for (size_t nu = igT[b - kolvoRegularNode]; nu < igT[b - kolvoRegularNode + 1]; nu++)
					{
						size_t pos_b = jgT[nu];
						if (a != pos_b)
						{
							if (pos_b > a)
								portrait[pos_b].insert(a);
							else
								portrait[a].insert(pos_b);
						}
					}
				}
				else
				{
					for (size_t mu = igT[a - kolvoRegularNode]; mu < igT[a - kolvoRegularNode + 1]; mu++)
					{
						size_t pos_a = jgT[mu];
						for (size_t nu = igT[b - kolvoRegularNode]; nu < igT[b - kolvoRegularNode + 1]; nu++)
						{
							size_t pos_b = jgT[nu];
							if (pos_b != pos_a)
							{
								if (pos_b > pos_a)
									portrait[pos_b].insert(pos_a);
								else
									portrait[pos_a].insert(pos_b);
							}
						}
					}
				}*/
			}
		}
	}

	size_t gg_size = 0;
	for (size_t i = 0; i < kolvoRegularNode; i++)
		gg_size += portrait[i].size();

	ig = new int[kolvoRegularNode + 1];
	di = new double[kolvoRegularNode];
	b = new double[kolvoRegularNode];
	q = new double[kolvoRegularNode];
	ig[0] = ig[1] = 0;
	for (size_t i = 0; i < kolvoRegularNode; i++)
	{
		di[i] = b[i] = q[i] = 0;
		//for (set<size_t>::iterator j = portrait[i].begin(); j != portrait[i].end(); ++j)
		//{
		//	slae.jg[tmp] = *j;
		//	tmp++;
		//}
		ig[i + 1] = ig[i] + portrait[i].size();
	}
	jg = new int[ig[kolvoRegularNode]];
	ggl = new double[ig[kolvoRegularNode]];
	ggu = new double[ig[kolvoRegularNode]];

	for (size_t i = 0; i < ig[kolvoRegularNode]; i++)
	{
		ggl[i] = 0;
		ggu[i] = 0;
	}

	size_t tmp = 0;
	for (size_t i = 0; i < kolvoRegularNode; i++)
	{
		for (set<size_t>::iterator j = portrait[i].begin(); j != portrait[i].end(); ++j)
		{
			jg[tmp] = *j;
			tmp++;
		}
		portrait[i].clear();
	}
	delete[] portrait;

	int i, j;
	ofstream igOut("ig.txt");
	ofstream jgOut("jg.txt");
	for (i = 0; i <= kolvoRegularNode; i++)
	{
		igOut << ig[i] << "   " << i << "   ";
		if (i != kolvoRegularNode)
			igOut << ig[i + 1] - ig[i] << endl;
	}
	cout << endl;
	for (j = 0; j < kolvoRegularNode; j++)
	{
		jgOut << "string " << j << ": ";
		for (i = ig[j]; i < ig[j + 1]; i++)
		{
			jgOut << jg[i] << " ";
		}
		jgOut << endl;
	}
}

void generatePortraitNesoglas()
{
	int kolvoRegularNode = xyz_points.size();
	int countLocalIndex = 8;
	set<size_t>* portrait = new set<size_t>[kolvoRegularNode];
	for (size_t k = 0; k < KE.size(); k++)
	{
		for (size_t i = 0; i < countLocalIndex; i++)
		{
			size_t a = KE[k].uzel[i];
			for (size_t j = 0; j < i; j++)
			{
				size_t b = KE[k].uzel[j];
				// Если оба узла не являются терминальными
				if (a < kolvoRegularNode && b < kolvoRegularNode)
				{
					if (b > a)
						portrait[b].insert(a);
					else
						portrait[a].insert(b);
				}
				else if (a >= kolvoRegularNode && b < kolvoRegularNode)
				{
					for (size_t mu = igT[a - kolvoRegularNode]; mu < igT[a - kolvoRegularNode + 1]; mu++)
					{
						size_t pos_a = jgT[mu];
						if (b != pos_a)
						{
							if (b > pos_a)
								portrait[b].insert(pos_a);
							else
								portrait[pos_a].insert(b);
						}
					}
				}
				else if (a < kolvoRegularNode && b >= kolvoRegularNode)
				{
					for (size_t nu = igT[b - kolvoRegularNode]; nu < igT[b - kolvoRegularNode + 1]; nu++)
					{
						size_t pos_b = jgT[nu];
						if (a != pos_b)
						{
							if (pos_b > a)
								portrait[pos_b].insert(a);
							else
								portrait[a].insert(pos_b);
						}
					}
				}
				else
				{
					for (size_t mu = igT[a - kolvoRegularNode]; mu < igT[a - kolvoRegularNode + 1]; mu++)
					{
						size_t pos_a = jgT[mu];
						for (size_t nu = igT[b - kolvoRegularNode]; nu < igT[b - kolvoRegularNode + 1]; nu++)
						{
							size_t pos_b = jgT[nu];
							if (pos_b != pos_a)
							{
								if (pos_b > pos_a)
									portrait[pos_b].insert(pos_a);
								else
									portrait[pos_a].insert(pos_b);
							}
						}
					}
				}
			}
		}
	}

	size_t gg_size = 0;
	for (size_t i = 0; i < kolvoRegularNode; i++)
		gg_size += portrait[i].size();

	ig = new int[kolvoRegularNode + 1];
	di = new double[kolvoRegularNode];
	b = new double[kolvoRegularNode];
	q = new double[kolvoRegularNode];
	ig[0] = ig[1] = 0;
	for (size_t i = 0; i < kolvoRegularNode; i++)
	{
		di[i] = b[i] = q[i] = 0;
		//for (set<size_t>::iterator j = portrait[i].begin(); j != portrait[i].end(); ++j)
		//{
		//	slae.jg[tmp] = *j;
		//	tmp++;
		//}
		ig[i + 1] = ig[i] + portrait[i].size();
	}
	jg = new int[ig[kolvoRegularNode]];
	ggl = new double[ig[kolvoRegularNode]];
	ggu = new double[ig[kolvoRegularNode]];

	for (size_t i = 0; i < ig[kolvoRegularNode]; i++)
	{
		ggl[i] = 0;
		ggu[i] = 0;
	}

	size_t tmp = 0;
	for (size_t i = 0; i < kolvoRegularNode; i++)
	{
		for (set<size_t>::iterator j = portrait[i].begin(); j != portrait[i].end(); ++j)
		{
			jg[tmp] = *j;
			tmp++;
		}
		portrait[i].clear();
	}
	delete[] portrait;

	int i, j;
	ofstream igOut("ig.txt");
	ofstream jgOut("jg.txt");
	for (i = 0; i <= kolvoRegularNode; i++)
	{
		igOut << ig[i] << "   " << i << "   ";
		if (i != kolvoRegularNode)
			igOut << ig[i + 1] - ig[i] << endl;
	}
	cout << endl;
	for (j = 0; j < kolvoRegularNode; j++)
	{
		jgOut << "string " << j << ": ";
		for (i = ig[j]; i < ig[j + 1]; i++)
		{
			jgOut << jg[i] << " ";
		}
		jgOut << endl;
	}
}

double analiticSolution(point goal)
{
	return 1 + goal.x + goal.y + goal.z;
	//return goal.x*goal.x + goal.y*goal.y + goal.z*goal.z;
	//return exp(goal.x + goal.y + goal.z);
}

double analiticSolution(double x, double y, double z)
{
	return analiticSolution(point(x, y, z));
}

double Lambda(int ielem)
{
	return 1;
}

double Gamma(int ielem)
{
	return 1;
}

double Func(int ielem, int node)
{
	return 1 + xyz_points[KE[ielem].uzel[node]].x + xyz_points[KE[ielem].uzel[node]].y + xyz_points[KE[ielem].uzel[node]].z;
	//return xyz_points[KE[ielem].uzel[node]].x*xyz_points[KE[ielem].uzel[node]].x + xyz_points[KE[ielem].uzel[node]].y*xyz_points[KE[ielem].uzel[node]].y +
	//	xyz_points[KE[ielem].uzel[node]].z*xyz_points[KE[ielem].uzel[node]].z - 6;
	//return analiticSolution(xyz_points[KE[ielem].uzel[node]].x, xyz_points[KE[ielem].uzel[node]].y, xyz_points[KE[ielem].uzel[node]].z) - 6;
	//return -2*analiticSolution(xyz_points[KE[ielem].uzel[node]].x, xyz_points[KE[ielem].uzel[node]].y, xyz_points[KE[ielem].uzel[node]].z);
}

void CreateLocalMatrixs(int ielem)
{
	double hXlocal, hYlocal, hZlocal;
	int i, j;
	hXlocal = xyz_points[KE[ielem].uzel[1]].x - xyz_points[KE[ielem].uzel[0]].x;
	hYlocal = xyz_points[KE[ielem].uzel[2]].y - xyz_points[KE[ielem].uzel[0]].y;
	hZlocal = xyz_points[KE[ielem].uzel[4]].z - xyz_points[KE[ielem].uzel[0]].z;
	double MWithoutGamma[8][8];
	for (i = 0; i < 8; i++)
	{
		short muI = (i % 2);
		short nuI = (i / 2) % 2;
		short etaI = (i / 4);
		localB[i] = 0;
		for (j = 0; j < 8; j++)
		{
			short muJ = (j % 2);
			short nuJ = (j / 2) % 2;
			short etaJ = (j / 4);

			localMatrix[i][j] = Lambda(ielem) * (hYlocal * hZlocal / hXlocal * G1[muI][muJ] * M1[nuI][nuJ] * M1[etaI][etaJ] +
				hXlocal * hZlocal / hYlocal * M1[muI][muJ] * G1[nuI][nuJ] * M1[etaI][etaJ] +
				hXlocal * hYlocal / hZlocal * M1[muI][muJ] * M1[nuI][nuJ] * G1[etaI][etaJ]);
			MWithoutGamma[i][j] = M1[muI][muJ] * M1[nuI][nuJ] * M1[etaI][etaJ] * hXlocal * hYlocal * hZlocal;
			localMatrix[i][j] += Gamma(ielem) * MWithoutGamma[i][j];

			localB[i] += MWithoutGamma[i][j] * Func(ielem, j);
		}
	}
}

void AddToMatrix(int posI, int posJ, double el)
{
	int tmp;
	if (posI == posJ)
	{
		di[posI] += el;
		return;
	}
	else
	{
		if (posI < posJ)
		{
			return;
		}
		for (tmp = ig[posI]; tmp < ig[posI + 1]; tmp++)
		{
			if (jg[tmp] == posJ)
			{
				ggl[tmp] += el;
				return;
			}
		}
	}
}

void Addition(int ielem)
{
	int i, j;
	int kolvoRegularNode = xyz_points.size();
	for (i = 0; i < 8; i++)
	{
		if (KE[ielem].uzel[i] < kolvoRegularNode)
		{
			b[KE[ielem].uzel[i]] += localB[i];
		}
		else
			for (int k = igT[KE[ielem].uzel[i] - kolvoRegularNode]; k < igT[KE[ielem].uzel[i] - kolvoRegularNode + 1]; k++)
			{
				b[jgT[k]] += localB[i] * ggT[k];
			}
		for (j = 0; j < 8; j++)
		{
			if (KE[ielem].uzel[i] < kolvoRegularNode)
			{
				if (KE[ielem].uzel[j] < kolvoRegularNode)
				{
					AddToMatrix(KE[ielem].uzel[i], KE[ielem].uzel[j], localMatrix[i][j]);
				}
				else
				{
					for (int k = igT[KE[ielem].uzel[j] - kolvoRegularNode]; k < igT[KE[ielem].uzel[j] - kolvoRegularNode + 1]; k++)
					{
						//posJ = jgT[k];
						//koefJ = ggT[k];
						AddToMatrix(KE[ielem].uzel[i], jgT[k], ggT[k] * localMatrix[i][j]);
					}
				}
			}
			else
			{
				if (KE[ielem].uzel[j] < kolvoRegularNode)
				{
					for (int k = igT[KE[ielem].uzel[i] - kolvoRegularNode]; k < igT[KE[ielem].uzel[i] - kolvoRegularNode + 1]; k++)
					{
						//posI = jgT[k];
						//koefI = ggT[k];
						AddToMatrix(jgT[k], KE[ielem].uzel[j], ggT[k] * localMatrix[i][j]);
					}
				}
				else
				{
					for (int k = igT[KE[ielem].uzel[i] - kolvoRegularNode]; k < igT[KE[ielem].uzel[i] - kolvoRegularNode + 1]; k++)
					{
						for (int l = igT[KE[ielem].uzel[j] - kolvoRegularNode]; l < igT[KE[ielem].uzel[j] - kolvoRegularNode + 1]; l++)
						{
							AddToMatrix(jgT[k], jgT[l], ggT[k] * ggT[l] * localMatrix[i][j]);
						}
					}
				}
			}
		}
	}
}

//intXorYorZ can be 0,1 or 2
void doEdge1(ofstream& outEdge1File, int intXorYorZ, int kolvoRegularNode, double* varNet1, int nVarNet1, double* varNet2, int nVarNet2, double uncnownValuePlosk)
{
	for (int iVar1 = 0; iVar1 < nVarNet1; iVar1 ++)
	{
		for (int iVar2 = 0; iVar2 < nVarNet2; iVar2++)
		{
			point goal;
			if (intXorYorZ == 0)
				goal = point(uncnownValuePlosk, varNet1[iVar1], varNet2[iVar2]);
			else
			{
				if (intXorYorZ == 1)
					goal = point(varNet1[iVar1], uncnownValuePlosk, varNet2[iVar2]);
				else
					goal = point(varNet1[iVar1], varNet2[iVar2], uncnownValuePlosk);
			}
			int k = indexXYZ(goal);
			di[k] = 1;
			b[k] = analiticSolution(goal);
			for (int m = ig[k]; m < ig[k + 1]; m++)
			{
				ggl[m] = 0;
			}
			for (int l = 0; l < kolvoRegularNode; l++)
			{
				for (int m = ig[l]; m < ig[l + 1]; m++)
				{
					if (k == jg[m])
					{
						ggu[m] = 0;
					}
				}
			}
			outEdge1File << k << '\t' << b[k] << endl;
		}
	}
}

void Edge1_not_sim(bool up, bool down, bool left, bool right, bool fore, bool behind)
{
	ofstream ku1("ku1.txt");
	int kolvoRegularNode = xyz_points.size();
	if (down)
		doEdge1(ku1, 2, kolvoRegularNode, xNet, nX, yNet, nY, zNet[0]);
	if (up)
		doEdge1(ku1, 2, kolvoRegularNode, xNet, nX, yNet, nY, zNet[nZ - 1]);
	if (left)
		doEdge1(ku1, 0, kolvoRegularNode, yNet, nY, zNet, nZ, xNet[0]);
	if (right)
		doEdge1(ku1, 0, kolvoRegularNode, yNet, nY, zNet, nZ, xNet[nX - 1]);
	if (fore)
		doEdge1(ku1, 1, kolvoRegularNode, xNet, nX, zNet, nZ, yNet[0]);
	if (behind)
		doEdge1(ku1, 1, kolvoRegularNode, xNet, nX, zNet, nZ, yNet[nY - 1]);
}

int findKE(int ind_nodes[4])
{
	for (size_t indKE = 0; indKE < KE.size(); indKE++)
	{
		int result[4] = {0,0,0,0};
		for (auto j = 0; j < 8; j++)
		{
			for (size_t k = 0; k < 4; k++)
			{
				if (KE[indKE].uzel[j] == ind_nodes[k])
				{
					result[k] = 1;
					break;
				}
			}
		}
		if (result[0] * result[1] * result[2] * result[3] == 1)
			return indKE;
	}

	cout << "Ошибка findKE. Не найден КЭ." << endl;
	system("pause");
	exit(1);
}

//intXorYorZ can be 0,1 or 2
void doEdge2(ofstream& outEdge2File, int intXorYorZ, int normalDirect, int kolvoRegularNode, double* varNet1, int nVarNet1, double* varNet2, int nVarNet2, double uncnownValuePlosk)
{
	double dUdn[4], h = 0.01, dh1, dh2;
	int indNodes[4];
	for (int iVar1 = 0; iVar1 < nVarNet1 - 1; iVar1++)
	{
		for (int iVar2 = 0; iVar2 < nVarNet2 - 1; iVar2++)
		{
			point goal;
			if (intXorYorZ == 0)
			{
				goal = point(uncnownValuePlosk, varNet1[iVar1], varNet2[iVar2]);
				dUdn[0] = normalDirect * (analiticSolution(uncnownValuePlosk + h, varNet1[iVar1], varNet2[iVar2]) - analiticSolution(uncnownValuePlosk - h, varNet1[iVar1], varNet2[iVar2])) / (2.0 * h);
				dUdn[1] = normalDirect * (analiticSolution(uncnownValuePlosk + h, varNet1[iVar1], varNet2[iVar2 + 1]) - analiticSolution(uncnownValuePlosk - h, varNet1[iVar1], varNet2[iVar2 + 1])) / (2.0 * h);
				dUdn[2] = normalDirect * (analiticSolution(uncnownValuePlosk + h, varNet1[iVar1 + 1], varNet2[iVar2]) - analiticSolution(uncnownValuePlosk - h, varNet1[iVar1 + 1], varNet2[iVar2])) / (2.0 * h);
				dUdn[3] = normalDirect * (analiticSolution(uncnownValuePlosk + h, varNet1[iVar1 + 1], varNet2[iVar2 + 1]) - analiticSolution(uncnownValuePlosk - h, varNet1[iVar1 + 1], varNet2[iVar2 + 1])) / (2.0 * h);
				indNodes[0] = indexXYZ(goal);
				indNodes[1] = indexXYZ(uncnownValuePlosk, varNet1[iVar1 + 1], varNet2[iVar2]);
				indNodes[2] = indexXYZ(uncnownValuePlosk, varNet1[iVar1], varNet2[iVar2 + 1]);
				indNodes[3] = indexXYZ(uncnownValuePlosk, varNet1[iVar1 + 1], varNet2[iVar2 + 1]);
			}
			else
			{
				if (intXorYorZ == 1)
				{
					goal = point(varNet1[iVar1], uncnownValuePlosk, varNet2[iVar2]);
					dUdn[0] = normalDirect * (analiticSolution(varNet1[iVar1], uncnownValuePlosk + h, varNet2[iVar2]) - analiticSolution(varNet1[iVar1], uncnownValuePlosk - h, varNet2[iVar2])) / (2.0 * h);
					dUdn[1] = normalDirect * (analiticSolution(varNet1[iVar1], uncnownValuePlosk + h, varNet2[iVar2 + 1]) - analiticSolution(varNet1[iVar1], uncnownValuePlosk - h, varNet2[iVar2 + 1])) / (2.0 * h);
					dUdn[2] = normalDirect * (analiticSolution(varNet1[iVar1 + 1], uncnownValuePlosk + h, varNet2[iVar2]) - analiticSolution(varNet1[iVar1 + 1], uncnownValuePlosk - h, varNet2[iVar2])) / (2.0 * h);
					dUdn[3] = normalDirect * (analiticSolution(varNet1[iVar1 + 1], uncnownValuePlosk + h, varNet2[iVar2 + 1]) - analiticSolution(varNet1[iVar1 + 1], uncnownValuePlosk - h, varNet2[iVar2 + 1])) / (2.0 * h);
					indNodes[0] = indexXYZ(goal);
					indNodes[1] = indexXYZ(varNet1[iVar1 + 1], uncnownValuePlosk, varNet2[iVar2]);
					indNodes[2] = indexXYZ(varNet1[iVar1], uncnownValuePlosk, varNet2[iVar2 + 1]);
					indNodes[3] = indexXYZ(varNet1[iVar1 + 1], uncnownValuePlosk, varNet2[iVar2 + 1]);
				}
				else
				{
					goal = point(varNet1[iVar1], varNet2[iVar2], uncnownValuePlosk);
					dUdn[0] = normalDirect * (analiticSolution(varNet1[iVar1], varNet2[iVar2], uncnownValuePlosk + h) - analiticSolution(varNet1[iVar1], varNet2[iVar2], uncnownValuePlosk - h)) / (2.0 * h);
					dUdn[1] = normalDirect * (analiticSolution(varNet1[iVar1], varNet2[iVar2 + 1], uncnownValuePlosk + h) - analiticSolution(varNet1[iVar1], varNet2[iVar2 + 1], uncnownValuePlosk - h)) / (2.0 * h);
					dUdn[2] = normalDirect * (analiticSolution(varNet1[iVar1 + 1], varNet2[iVar2], uncnownValuePlosk + h) - analiticSolution(varNet1[iVar1 + 1], varNet2[iVar2], uncnownValuePlosk - h)) / (2.0 * h);
					dUdn[3] = normalDirect * (analiticSolution(varNet1[iVar1 + 1], varNet2[iVar2 + 1], uncnownValuePlosk + h) - analiticSolution(varNet1[iVar1 + 1], varNet2[iVar2 + 1], uncnownValuePlosk - h)) / (2.0 * h);
					indNodes[0] = indexXYZ(goal);
					indNodes[1] = indexXYZ(varNet1[iVar1 + 1], varNet2[iVar2], uncnownValuePlosk);
					indNodes[2] = indexXYZ(varNet1[iVar1], varNet2[iVar2 + 1], uncnownValuePlosk);
					indNodes[3] = indexXYZ(varNet1[iVar1 + 1], varNet2[iVar2 + 1], uncnownValuePlosk);
				}
			}
			int indKE = findKE(indNodes);
			dh1 = varNet1[iVar1 + 1] - varNet1[iVar1];
			dh2 = varNet2[iVar2 + 1] - varNet2[iVar2];

			for (size_t linux = 0; linux < 4; linux++)
			{
				double value = 0;
				for (size_t linux2 = 0; linux2 < 4; linux2 ++)
				{
					value += dUdn[linux2] * M2[linux][linux2] * dh1 * dh2 / 36.0 * Lambda(indKE);
				}
				b[indNodes[linux]] += value;
				outEdge2File << setw(10) << indNodes[linux] << setw(15) << value << endl;
			}
		}
	}
}

void Edge2_not_sim(bool up, bool down, bool left, bool right, bool fore, bool behind)
{
	ofstream ku2("ku2.txt");
	int kolvoRegularNode = xyz_points.size();
	if (up)
		doEdge2(ku2, 2, 1, kolvoRegularNode, xNet, nX, yNet, nY, zNet[nZ - 1]);
	if (down)
		doEdge2(ku2, 2, -1, kolvoRegularNode, xNet, nX, yNet, nY, zNet[0]);
	if (left)
		doEdge2(ku2, 0, -1, kolvoRegularNode, yNet, nY, zNet, nZ, xNet[0]);
	if (right)
		doEdge2(ku2, 0, 1, kolvoRegularNode, yNet, nY, zNet, nZ, xNet[nX - 1]);
	if (fore)
		doEdge2(ku2, 1, -1, kolvoRegularNode, xNet, nX, zNet, nZ, yNet[0]);
	if (behind)
		doEdge2(ku2, 1, 1, kolvoRegularNode, xNet, nX, zNet, nZ, yNet[nY - 1]);
}

//intXorYorZ can be 0,1 or 2
void doEdge3(ofstream& outEdge3File, int intXorYorZ, int normalDirect, int kolvoRegularNode, double* varNet1, int nVarNet1, double* varNet2, int nVarNet2, double uncnownValuePlosk)
{
	double dUdn[4], h = 0.01, dh1, dh2, ubetta[4], betta;
	int indNodes[4];
	for (int iVar1 = 0; iVar1 < nVarNet1 - 1; iVar1++)
	{
		for (int iVar2 = 0; iVar2 < nVarNet2 - 1; iVar2++)
		{
			point goal;
			if (intXorYorZ == 0)
			{
				goal = point(uncnownValuePlosk, varNet1[iVar1], varNet2[iVar2]);
				dUdn[0] = normalDirect * (analiticSolution(uncnownValuePlosk + h, varNet1[iVar1], varNet2[iVar2]) - analiticSolution(uncnownValuePlosk - h, varNet1[iVar1], varNet2[iVar2])) / (2.0 * h);
				dUdn[1] = normalDirect * (analiticSolution(uncnownValuePlosk + h, varNet1[iVar1], varNet2[iVar2 + 1]) - analiticSolution(uncnownValuePlosk - h, varNet1[iVar1], varNet2[iVar2 + 1])) / (2.0 * h);
				dUdn[2] = normalDirect * (analiticSolution(uncnownValuePlosk + h, varNet1[iVar1 + 1], varNet2[iVar2]) - analiticSolution(uncnownValuePlosk - h, varNet1[iVar1 + 1], varNet2[iVar2])) / (2.0 * h);
				dUdn[3] = normalDirect * (analiticSolution(uncnownValuePlosk + h, varNet1[iVar1 + 1], varNet2[iVar2 + 1]) - analiticSolution(uncnownValuePlosk - h, varNet1[iVar1 + 1], varNet2[iVar2 + 1])) / (2.0 * h);
				indNodes[0] = indexXYZ(goal);
				indNodes[1] = indexXYZ(uncnownValuePlosk, varNet1[iVar1 + 1], varNet2[iVar2]);
				indNodes[2] = indexXYZ(uncnownValuePlosk, varNet1[iVar1], varNet2[iVar2 + 1]);
				indNodes[3] = indexXYZ(uncnownValuePlosk, varNet1[iVar1 + 1], varNet2[iVar2 + 1]);
			}
			else
			{
				if (intXorYorZ == 1)
				{
					goal = point(varNet1[iVar1], uncnownValuePlosk, varNet2[iVar2]);
					dUdn[0] = normalDirect * (analiticSolution(varNet1[iVar1], uncnownValuePlosk + h, varNet2[iVar2]) - analiticSolution(varNet1[iVar1], uncnownValuePlosk - h, varNet2[iVar2])) / (2.0 * h);
					dUdn[1] = normalDirect * (analiticSolution(varNet1[iVar1], uncnownValuePlosk + h, varNet2[iVar2 + 1]) - analiticSolution(varNet1[iVar1], uncnownValuePlosk - h, varNet2[iVar2 + 1])) / (2.0 * h);
					dUdn[2] = normalDirect * (analiticSolution(varNet1[iVar1 + 1], uncnownValuePlosk + h, varNet2[iVar2]) - analiticSolution(varNet1[iVar1 + 1], uncnownValuePlosk - h, varNet2[iVar2])) / (2.0 * h);
					dUdn[3] = normalDirect * (analiticSolution(varNet1[iVar1 + 1], uncnownValuePlosk + h, varNet2[iVar2 + 1]) - analiticSolution(varNet1[iVar1 + 1], uncnownValuePlosk - h, varNet2[iVar2 + 1])) / (2.0 * h);
					indNodes[0] = indexXYZ(goal);
					indNodes[1] = indexXYZ(varNet1[iVar1 + 1], uncnownValuePlosk, varNet2[iVar2]);
					indNodes[2] = indexXYZ(varNet1[iVar1], uncnownValuePlosk, varNet2[iVar2 + 1]);
					indNodes[3] = indexXYZ(varNet1[iVar1 + 1], uncnownValuePlosk, varNet2[iVar2 + 1]);
				}
				else
				{
					goal = point(varNet1[iVar1], varNet2[iVar2], uncnownValuePlosk);
					dUdn[0] = normalDirect * (analiticSolution(varNet1[iVar1], varNet2[iVar2], uncnownValuePlosk + h) - analiticSolution(varNet1[iVar1], varNet2[iVar2], uncnownValuePlosk - h)) / (2.0 * h);
					dUdn[1] = normalDirect * (analiticSolution(varNet1[iVar1], varNet2[iVar2 + 1], uncnownValuePlosk + h) - analiticSolution(varNet1[iVar1], varNet2[iVar2 + 1], uncnownValuePlosk - h)) / (2.0 * h);
					dUdn[2] = normalDirect * (analiticSolution(varNet1[iVar1 + 1], varNet2[iVar2], uncnownValuePlosk + h) - analiticSolution(varNet1[iVar1 + 1], varNet2[iVar2], uncnownValuePlosk - h)) / (2.0 * h);
					dUdn[3] = normalDirect * (analiticSolution(varNet1[iVar1 + 1], varNet2[iVar2 + 1], uncnownValuePlosk + h) - analiticSolution(varNet1[iVar1 + 1], varNet2[iVar2 + 1], uncnownValuePlosk - h)) / (2.0 * h);
					indNodes[0] = indexXYZ(goal);
					indNodes[1] = indexXYZ(varNet1[iVar1 + 1], varNet2[iVar2], uncnownValuePlosk);
					indNodes[2] = indexXYZ(varNet1[iVar1], varNet2[iVar2 + 1], uncnownValuePlosk);
					indNodes[3] = indexXYZ(varNet1[iVar1 + 1], varNet2[iVar2 + 1], uncnownValuePlosk);
				}
			}

			int indKE = findKE(indNodes);

			dh1 = varNet1[iVar1 + 1] - varNet1[iVar1];
			dh2 = varNet2[iVar2 + 1] - varNet2[iVar2];
			betta = -dUdn[0] * Lambda(indKE);

			for (size_t i = 0; i < 4; i++)
			{
				ubetta[i] = analiticSolution(xyz_points[indNodes[i]]) - 1;
			}

			for (size_t linux = 0; linux < 4; linux++)
			{
				double value = 0;
				for (size_t linux2 = 0; linux2 < 4; linux2++)
				{
					value += ubetta[linux2] * M2[linux][linux2] * dh1 * dh2 / 36.0;
					AddToMatrix(indNodes[linux], indNodes[linux2], betta * M2[linux][linux2] * dh1 * dh2 / 36.0);
				}
				b[indNodes[linux]] += betta * value;
				outEdge3File << setw(10) << indNodes[linux] << setw(15) << value << endl;
			}
		}
	}
}

void Edge3_not_sim(bool up, bool down, bool left, bool right, bool fore, bool behind)
{
	ofstream ku3("ku3.txt");
	int kolvoRegularNode = xyz_points.size();
	if (up)
		doEdge3(ku3, 2, 1, kolvoRegularNode, xNet, nX, yNet, nY, zNet[nZ - 1]);
	if (down)
		doEdge3(ku3, 2, -1, kolvoRegularNode, xNet, nX, yNet, nY, zNet[0]);
	if (left)
		doEdge3(ku3, 0, -1, kolvoRegularNode, yNet, nY, zNet, nZ, xNet[0]);
	if (right)
		doEdge3(ku3, 0, 1, kolvoRegularNode, yNet, nY, zNet, nZ, xNet[nX - 1]);
	if (fore)
		doEdge3(ku3, 1, -1, kolvoRegularNode, xNet, nX, zNet, nZ, yNet[0]);
	if (behind)
		doEdge3(ku3, 1, 1, kolvoRegularNode, xNet, nX, zNet, nZ, yNet[nY - 1]);
}

void GenerateMatrix()
{
	int ielem, i;
	int kolvoRegularNode = xyz_points.size();
	for (ielem = 0; ielem < KE.size(); ielem++)
	{
		CreateLocalMatrixs(ielem);
		Addition(ielem);
	}

	//Edge2_not_sim(0, 0, 1, 1, 0, 0);
	//Edge3_not_sim(0, 0, 0, 0, 1, 1);
	for (i = 0; i < ig[kolvoRegularNode]; i++)
	{
		ggu[i] = ggl[i];
	}
	//Edge1_not_sim(1, 1, 0, 0, 0, 0);
	Edge1_not_sim(1, 1, 1, 1, 1, 1);
}

void mult(double* res, double* v)
{
	int kolvoRegularNode = xyz_points.size();
	for (int i = 0; i < kolvoRegularNode; i++)
		res[i] = 0;
	for (int i = 0; i < kolvoRegularNode; i++)
	{
		for (int j = ig[i]; j < ig[i + 1]; j++)
		{
			res[i] += ggl[j] * v[jg[j]];
			res[jg[j]] += ggu[j] * v[i];
		}
		res[i] += di[i] * v[i];
	}
}

double ScalarMult(double* v1, double* v2)
{
	int i;
	double result;
	int kolvoRegularNode = xyz_points.size();
	result = 0;
	for (i = 0; i < kolvoRegularNode; i++)
	{
		result += v1[i] * v2[i];
	}
	return result;
}

//	Умножение матрицы на вектор
void MultMatrixOnVector(double* in, double* out)
{
	int i, j;
	double* out1;
	int kolvoRegularNode = xyz_points.size();
	out1 = new double[kolvoRegularNode];
	for (i = 0; i < kolvoRegularNode; i++)
	{
		out1[i] = di[i] * in[i];
		for (j = ig[i]; j < ig[i + 1]; j++)
		{
			out1[i] += ggl[j] * in[jg[j]];
			out1[jg[j]] += ggu[j] * in[i];
		}
	}
	for (i = 0; i < kolvoRegularNode; i++)
		out[i] = out1[i];
	delete[] out1;
}

void calcPogreshnost(ofstream& output)
{
	int i;
	double sumCh = 0, sumZn = 0;
	output << setw(10) << "x" << setw(10) << "y" << setw(10) << "z" << setw(18) << "analitic" << setw(18) << "solution" << setw(18) << "pogreshn" << endl;
	for (i = 0; i < xyz_points.size(); i++)
	{
		sumCh += (q[i] - analiticSolution(xyz_points[i])) * (q[i] - analiticSolution(xyz_points[i]));
		output << setw(10) << xyz_points[i].x << setw(10) << xyz_points[i].y << setw(10) << xyz_points[i].z << setw(18) << analiticSolution(xyz_points[i]) << setw(18) << q[i] << setw(18) << q[i] - analiticSolution(xyz_points[i]) << endl;
		sumZn += analiticSolution(xyz_points[i]) * analiticSolution(xyz_points[i]);
	}
	output << "KE number = " << KE.size() << endl;
	output << "Nodes number = " << xyz_points.size() << endl;
	output << endl << "Otnositelnaia pogreshnost = " << sqrt(sumCh / sumZn);
}

void runLOS()
{
	int maxiter = 10000, i;
	double alfa, alfachisl, alfaznam, beta, betachisl, betaznam, checkE, epsMSG = 1e-16;

	int kolvoRegularNode = xyz_points.size();
	double* r = new double[kolvoRegularNode];
	double* s = new double[kolvoRegularNode];
	double* z = new double[kolvoRegularNode];
	double* p = new double[kolvoRegularNode];
	double* rout = new double[kolvoRegularNode];

	for (i = 0; i < kolvoRegularNode; i++)
	{
		s[i] = rout[i] = r[i] = q[i] = z[i] = p[i] = 0;
	}
	MultMatrixOnVector(q, r);
	for (i = 0; i < kolvoRegularNode; i++)
	{
		r[i] = b[i] - r[i];
		z[i] = r[i];
	}
	MultMatrixOnVector(z, p);
	checkE = sqrt(ScalarMult(r, r) / ScalarMult(b, b));
	//startNeviazka = checkE = ScalarMult(r, r);
	for (int iter = 0; iter < maxiter && checkE >= epsMSG; iter++)
	{
		alfachisl = ScalarMult(p, r);
		alfaznam = ScalarMult(p, p);
		alfa = alfachisl / alfaznam;
		for (i = 0; i < kolvoRegularNode; i++)
		{
			q[i] = q[i] + alfa * z[i];
			r[i] = r[i] - alfa * p[i];
		}
		MultMatrixOnVector(r, rout);
		betachisl = ScalarMult(p, rout);
		betaznam = ScalarMult(p, p);
		beta = -betachisl / betaznam;
		for (i = 0; i < kolvoRegularNode; i++)
		{
			z[i] = r[i] + beta * z[i];
			p[i] = rout[i] + beta * p[i];
		}
		checkE = sqrt(ScalarMult(r, r) / ScalarMult(b, b));
	}
}

void outputWithoutOptimisation()
{
	ofstream fileXY("xyz.txt");
	ofstream fileNvtr("nvtr.txt");

	for (int k = 0; k < nZ; k++)
	{
		for (int j = 0; j < nY; j++)
		{
			for (int i = 0; i < nX; i++)
			{
				xyz_points.push_back(point(xNet[i], yNet[j], zNet[k]));
				fileXY << xNet[i] << " " << yNet[j] << " " << zNet[k] << endl;
			}
		}
	}

	for (int zI = 0; zI < nZ - 1; ++zI)
	{
		for (int yI = 0; yI < nY - 1; ++yI)
		{
			for (int xI = 0; xI < nX - 1; ++xI)
			{
				nvtr tempNvtr;
				tempNvtr.uzel[0] = indexXYZ(xNet[xI], yNet[yI], zNet[zI]);
				tempNvtr.uzel[1] = indexXYZ(xNet[xI + 1], yNet[yI], zNet[zI]);
				tempNvtr.uzel[2] = indexXYZ(xNet[xI], yNet[yI + 1], zNet[zI]);
				tempNvtr.uzel[3] = indexXYZ(xNet[xI + 1], yNet[yI + 1], zNet[zI]);

				tempNvtr.uzel[4] = indexXYZ(xNet[xI], yNet[yI], zNet[zI + 1]);
				tempNvtr.uzel[5] = indexXYZ(xNet[xI + 1], yNet[yI], zNet[zI + 1]);
				tempNvtr.uzel[6] = indexXYZ(xNet[xI], yNet[yI + 1], zNet[zI + 1]);
				tempNvtr.uzel[7] = indexXYZ(xNet[xI + 1], yNet[yI + 1], zNet[zI + 1]);

				tempNvtr.numberField = 0;

				KE.push_back(tempNvtr);
				for (int i = 0; i < 8; i++)
				{
					fileNvtr << tempNvtr.uzel[i] << " ";
				}
				fileNvtr << tempNvtr.numberField << endl;
			}
		}
	}
}

int findAreaNumber(int nodes[])
{
	int i;
	for (i = 0; i < sreda.size(); i++)
	{
		if (xyz_points[nodes[0]].x >= sreda[i].x1 && xyz_points[nodes[1]].x <= sreda[i].x2 &&
			xyz_points[nodes[0]].y >= sreda[i].y1 && xyz_points[nodes[2]].y <= sreda[i].y2 &&
			xyz_points[nodes[0]].z >= sreda[i].z1 && xyz_points[nodes[4]].z <= sreda[i].z2)
			return i;
	}
	cout << "Ошибка в FindAreaNumber: не найдена подобласть." << endl;
	system("pause");
	exit(1);
	return -1;
}

void outputKEandXYZ()
{
	ofstream fileXY("xyz.txt");
	ofstream fileNvtr("nvtr.txt");
	ofstream viewNet("viewNet.txt");
	vector<point> ncPoint;

	//формируем массивы регулярных и терминальных вершин
	for (int k = 0; k < nZ; k++)
		for (int j = 0; j < nY; j++)
			for (int i = 0; i < nX; i++)
			{
				if (matrixNode[i][j][k] != 'Y')
					if ((k == 0 && matrixNode[i][j][k] == 'F' || k == nZ - 1 && matrixNode[i][j][k] == 'B'
							|| i == 0 && matrixNode[i][j][k] == 'A' || i == nX - 1 && matrixNode[i][j][k] == 'D'
							|| j == 0 && matrixNode[i][j][k] == 'S' || j == nY - 1 && matrixNode[i][j][k] == 'W')
						|| matrixNode[i][j][k] == 'V' || matrixNode[i][j][k] == 'R')
					{
						xyz_points.push_back(point(xNet[i], yNet[j], zNet[k]));
					}
					else
					{
						ncPoint.push_back(point(xNet[i], yNet[j], zNet[k]));
					}
			}

	//добавляем в КОНЕЦ массива вершин терминальные вершины
	for (vector<point>::iterator it = ncPoint.begin(); it < ncPoint.end(); it++)
	{
		xyz_points.push_back(*it);
	}

	//количество столбцов в T
	nColT = ncPoint.size();
	kolvoRegularNode = xyz_points.size() - nColT;

	//формируем файл xy.txt
	fileXY << xyz_points.size() << " " << kolvoRegularNode << endl;
	for (int i = 0; i < xyz_points.size(); i++)
	{
		fileXY << xyz_points[i].x << " " << xyz_points[i].y << " " << xyz_points[i].z << endl;
	}

	for (int iZ = 0; iZ < nZ; iZ++)
	{
		for (int iY = nY - 1; iY >= 0; iY--)
		{
			for (int iX = 0; iX < nX; iX++)
			{
				viewNet << matrixNode[iX][iY][iZ] << " ";
			}
			viewNet << endl;
		}
		viewNet << "-----------------------" << endl;
	}
	viewNet.close();

	//формируем структуру КЭ
	nvtr tempNvtr;
	for (int iZ = 0; iZ < nZ - 1; iZ++)
		for (int iY = 0; iY < nY - 1; iY++)
			for (int iX = 0; iX < nX - 1; iX++)
			{
				if (matrixNode[iX][iY][iZ] == 'V' || matrixNode[iX][iY][iZ] == 'S' || matrixNode[iX][iY][iZ] == 'A'
					|| matrixNode[iX][iY][iZ] == 'R' || matrixNode[iX][iY][iZ] == 'F')
				{
					tempNvtr.uzel[0] = indexXYZ(xNet[iX], yNet[iY], zNet[iZ]);
					for (int k = iX + 1; k < nX; k++)
					{
						if (matrixNode[k][iY][iZ] == 'S' || matrixNode[k][iY][iZ] == 'R' || matrixNode[k][iY][iZ] == 'D'
							|| matrixNode[k][iY][iZ] == 'V' || matrixNode[k][iY][iZ] == 'F')
						{
							tempNvtr.uzel[1] = indexXYZ(xNet[k], yNet[iY], zNet[iZ]);
							for (int t = iY + 1; t < nY; t++)
							{
								if (matrixNode[k][t][iZ] == 'W' || matrixNode[k][t][iZ] == 'R' || matrixNode[k][t][iZ] == 'D'
									|| matrixNode[k][t][iZ] == 'V' || matrixNode[k][t][iZ] == 'F')
								{
									tempNvtr.uzel[3] = indexXYZ(xNet[k], yNet[t], zNet[iZ]);
									tempNvtr.uzel[2] = indexXYZ(xNet[iX], yNet[t], zNet[iZ]);

									for (int h = iZ + 1; h < nZ; h++)
									{
										if (matrixNode[iX][iY][h] == 'S' || matrixNode[iX][iY][h] == 'R' || matrixNode[iX][iY][h] == 'A'
											|| matrixNode[iX][iY][h] == 'V' || matrixNode[iX][iY][h] == 'B')
										{
											tempNvtr.uzel[4] = indexXYZ(xNet[iX], yNet[iY], zNet[h]);
											tempNvtr.uzel[5] = indexXYZ(xNet[k], yNet[iY], zNet[h]);
											tempNvtr.uzel[6] = indexXYZ(xNet[iX], yNet[t], zNet[h]);
											tempNvtr.uzel[7] = indexXYZ(xNet[k], yNet[t], zNet[h]);
											break;
										}
									}
									break;
								}
							}
							break;
						}
					}
					tempNvtr.numberField = findAreaNumber(tempNvtr.uzel);
					KE.push_back(tempNvtr);
				}
			}

	//формируем файл nvtr.txt
	fileNvtr << KE.size() << endl;
	for (int i = 0; i < KE.size(); i++)
	{
		for (int j = 0; j < 8; j++)
		{
			fileNvtr << KE[i].uzel[j] << " ";
		}
		fileNvtr << KE[i].numberField << endl;
	}
}

locateOfPoint FindLocate(point sample)
{
	locateOfPoint a;
	a.i = -1;
	a.j = -1;
	a.k = -1;
	for (int i = 0; i < nX; i++)
	{
		if (xNet[i] == sample.x)
			a.i = i;
	}
	for (int j = 0; j < nY; j++)
	{
		if (yNet[j] == sample.y)
			a.j = j;
	}
	for (int k = 0; k < nY; k++)
	{
		if (zNet[k] == sample.z)
			a.k = k;
	}
	if (a.i == -1 || a.j == -1 || a.k == -1)
	{
		cout << "Ошибка FindLocate: не найдена точка" << endl;
		exit(1);
	}
	return a;
}

void sigmTChain(int nTermNode, int startOfChain, double mnojT, set<int> &visitedNodes)
{
	int i, j;
	for each (neighbor var in tmpSigm[startOfChain].neighbors)
	{
		//if (!WithoutPerehlest)
			//for (i = 0; i < sigmT[nTermNode].neighbors.size(); i++)
			//{
			//	if (sigmT[nTermNode].neighbors[i] == sigmT[startOfChain].neighbors[j])
			//	{
			//		cout << "Обнаружена цепочка терминальных узлов в " << nTermNode << endl;
			//		//system("pause");
			//		//exit(1);
			//	}
			//}
		if (visitedNodes.find(var.index) != visitedNodes.end()) {
			logger << "found circle " << var.index << " in node " << startOfChain << endl;
			continue;
		}
		visitedNodes.insert(var.index);
		if (var.index >= kolvoRegularNode)
		{
			sigmTChain(nTermNode, var.index - kolvoRegularNode, mnojT * var.weight, visitedNodes);
		} else {
			sigmNewT[nTermNode].neighbors.insert(neighbor(var.index, mnojT*var.weight));
		}
	}
}

void genT3D() {
	ofstream outT("Tmatrix.txt");
	igT = new int[nColT + 1];

	//формируем вспомогательную сигм-структуру
	//sigmNewT = new sigmStruct3D[nColT];
	sigmNewT.resize(nColT);
	tmpSigm = new sigmStruct3D[nColT];

	for (int i = kolvoRegularNode; i < xyz_points.size(); i++)
	{
		int inc = i - kolvoRegularNode;
		tmpSigm[inc].terminalNode = i;
		locateOfPoint term = FindLocate(xyz_points[i]);
		if (newNodes[term.i][term.j][term.k].test(IS_REGULAR))
			throw new exception("wrong terminal node in genT3D");
		if (hasXLines(newNodes[term.i][term.j][term.k]) && hasYLines(newNodes[term.i][term.j][term.k])) {
			bool downResult = hasDown(newNodes[term.i][term.j][term.k]);
			bool upResult = hasUp(newNodes[term.i][term.j][term.k]);
			if (DEBUG && (downResult && upResult || !downResult && !upResult)) {
				logError("something goes wrong with term nodes");
			}
			if (upResult) {

			}
		}
		if (hasXLines(newNodes[term.i][term.j][term.k]))
		{
			int neibLeft;
			for (neibLeft = term.i - 1; neibLeft >= 0; neibLeft--)
			{
				if (!newNodes[neibLeft][term.j][term.k].none()) break;
			}
			int neibRight;
			for (neibRight = term.i + 1; neibRight < nX; neibRight++)
			{
				if (!newNodes[neibRight][term.j][term.k].none()) break;
			}
			double length = xNet[neibRight] - xNet[neibLeft];
			int iL = indexXYZ(xNet[neibLeft], yNet[term.j], zNet[term.k]);
			int iR = indexXYZ(xNet[neibRight], yNet[term.j], zNet[term.k]);
			tmpSigm[inc].neighbors.insert(neighbor(iL, (xNet[neibRight] - xNet[term.i]) / length));
			tmpSigm[inc].neighbors.insert(neighbor(iR, (xNet[term.i] - xNet[neibLeft]) / length));
		}
		if (hasYLines(newNodes[term.i][term.j][term.k]))
		{
			int neibFore;
			for (neibFore = term.j + 1; neibFore < nY; neibFore++)
			{
				if (!newNodes[term.i][neibFore][term.k].none()) break;
			}
			int neibBack;
			for (neibBack = term.j - 1; neibBack >= 0; neibBack--)
			{
				if (!newNodes[term.i][neibBack][term.k].none()) break;
			}
			double length = yNet[neibBack] - yNet[neibFore];
			int iF = indexXYZ(xNet[term.i], yNet[neibFore], zNet[term.k]);
			int iB = indexXYZ(xNet[term.i], yNet[neibBack], zNet[term.k]);
			tmpSigm[inc].neighbors.insert(neighbor(iF, (yNet[neibBack] - yNet[term.j]) / length));
			tmpSigm[inc].neighbors.insert(neighbor(iB, (yNet[term.j] - yNet[neibFore]) / length));
		}
		if (hasZLines(newNodes[term.i][term.j][term.k]))
		{
			int neibDown;
			for (neibDown = term.k + 1; neibDown < nZ; neibDown++)
			{
				if (!newNodes[term.i][term.j][neibDown].none()) break;
			}
			int neibUp;
			for (neibUp = term.k - 1; neibUp >= 0; neibUp--)
			{
				if (!newNodes[term.i][term.j][neibUp].none()) break;
			}
			double length = yNet[neibUp] - yNet[neibDown];
			int iD = indexXYZ(xNet[term.i], yNet[term.j], zNet[neibDown]);
			int iU = indexXYZ(xNet[term.i], yNet[term.j], zNet[neibUp]);
			tmpSigm[inc].neighbors.insert(neighbor(iD, (zNet[neibUp] - zNet[term.k]) / length));
			tmpSigm[inc].neighbors.insert(neighbor(iU, (zNet[term.k] - zNet[neibDown]) / length));
		}
		if (tmpSigm[inc].neighbors.empty())
			throw new exception("cannot find neighbors for tmpSigm");
	}

	for (int i = 0; i < nColT; i++)
	{
		set<int> visitedNodes;
		visitedNodes.insert(tmpSigm[i].terminalNode);
		sigmNewT[i].terminalNode = tmpSigm[i].terminalNode;
		logger << "Start chaining " << tmpSigm[i].terminalNode << " node" << endl;
		for each (neighbor var in tmpSigm[i].neighbors)
		{
			visitedNodes.insert(var.index);
			if (var.index >= kolvoRegularNode) {
				sigmTChain(i, var.index - kolvoRegularNode, var.weight, visitedNodes);
			}
			else {
				sigmNewT[i].neighbors.insert(var);
			}
		}
		//check
		double sum = 0;
		for each (neighbor var in sigmNewT[i].neighbors)
		{
			sum += var.weight;
		}
		if (!MkaUtils::equals(sum, 1)) {
			logError("Incorrect chain of matrix T, weights sum not equals to 1");
		}
	}

	for (int i = 0; i < nColT; i++)
	{
		outT << sigmNewT[i].terminalNode << " ";
		for each (neighbor var in sigmNewT[i].neighbors)
		{
			outT << var.index << " " << var.weight << " ";
		}
		outT << endl;
	}
	igT[0] = 0;
	for (int i = 0; i < nColT; i++)
	{
		igT[i + 1] = igT[i] + sigmNewT[i].neighbors.size();
	}
	ggT = new double[igT[nColT]];
	jgT = new int[igT[nColT]];
	for (int i = 0, j = 0; i < nColT && j < igT[nColT]; i++)
	{
		for each (neighbor var in sigmNewT[i].neighbors)
		{
			jgT[j] = var.index;
			ggT[j] = var.weight;
			j++;
		}
	}

	outT << endl << "T" << endl;
	for (int i = 0; i <= nColT; i++)
	{
		outT << igT[i] << " ";
	}
	outT << endl;
	for (int i = 0; i < igT[nColT]; i++)
	{
		outT << jgT[i] << " ";
	}
	outT << endl;
	for (int i = 0; i < igT[nColT]; i++)
	{
		outT << ggT[i] << " ";
	}
	outT << endl;
}

void Output()
{
	int i, j, t, k;
	ofstream fileXY("xyz.txt");
	ofstream fileNvtr("nvtr.txt");
	vector<point> ncPoint;

	//формируем массивы регулярных и терминальных вершин
	for (t = 0; t < nZ; t++)
	{
		for (j = 0; j < nY; j++)
		{
			for (i = 0; i < nX; i++)
			{
				if (DEBUG && !newNodes[i][j][t].test(IS_REGULAR) && hasAll(newNodes[i][j][t]))
					throw new exception("it is happened!");
				if (newNodes[i][j][t].test(IS_REGULAR) || hasAll(newNodes[i][j][t]))	//может быть и терминальным, добавить проверку
				{
					xyz_points.push_back(point(xNet[i], yNet[j],zNet[t]));
				}
				else if (!newNodes[i][j][t].none())
				{
					ncPoint.push_back(point(xNet[i], yNet[j], zNet[t]));
				}
			}
		}
	}

	//добавляем в конец массива вершин терминальные вершины
	for (vector<point>::iterator it = ncPoint.begin(); it < ncPoint.end(); it++)
	{
		xyz_points.push_back(*it);
	}

	//количество столбцов в T
	nColT = ncPoint.size();
	kolvoRegularNode = xyz_points.size() - nColT;

	//формируем файл xyz.txt
	fileXY << xyz_points.size() << " " << kolvoRegularNode << endl;
	for (i = 0; i < xyz_points.size(); i++)
	{
		fileXY << xyz_points[i].x << " " << xyz_points[i].y << " " << xyz_points[i].z << endl;
	}

	//формируем структуру КЭ
	nvtr tempNvtr;
	for (int z = 0; z < nZ - 1; z++)
	{
		for (j = 0; j < nY - 1; j++)
		{
			for (i = 0; i < nX - 1; i++)
			{
				if (canOptimize(1, 1, 1, newNodes[i][j][z]))
				{
					int nextX = findNextX(1, 1, 1, i, nX, j, z);
					int nextY = findNextY(1, 1, 1, j, nY, i, z);
					int nextZ = findNextZ(1, 1, 1, z, nZ, i, j);

					tempNvtr.uzel[0] = indexXYZ(xNet[i], yNet[j], zNet[z]);
					tempNvtr.uzel[1] = indexXYZ(xNet[nextX], yNet[j], zNet[z]);
					tempNvtr.uzel[2] = indexXYZ(xNet[i], yNet[nextY], zNet[z]);
					tempNvtr.uzel[3] = indexXYZ(xNet[nextX], yNet[nextY], zNet[z]);

					tempNvtr.uzel[4] = indexXYZ(xNet[i], yNet[j], zNet[nextZ]);
					tempNvtr.uzel[5] = indexXYZ(xNet[nextX], yNet[j], zNet[nextZ]);
					tempNvtr.uzel[6] = indexXYZ(xNet[i], yNet[nextY], zNet[nextZ]);
					tempNvtr.uzel[7] = indexXYZ(xNet[nextX], yNet[nextY], zNet[nextZ]);

					tempNvtr.numberField = 1; //FindAreaNumber(tempNvtr.uzel);
					KE.push_back(tempNvtr);
				}
			}
		}
	}

	//формируем файл nvtr.txt
	fileNvtr << KE.size()
		<< "\t" << xNet[0] << "\t" << xNet[nX - 1]
		<< "\t" << yNet[0] << "\t" << yNet[nY - 1]
		<< "\t" << zNet[0] << "\t" << zNet[nZ - 1] << endl;
	for (i = 0; i < KE.size(); i++)
	{
		for (j = 0; j < 8; j++)
		{
			fileNvtr << KE[i].uzel[j] << " ";
		}
		fileNvtr << endl;
	}
}

int main(int argc, char* argv[])
{
	setlocale(LC_ALL, "rus");
	inputConfig();
	inputNet();

	initNet(xNet, nX, yNet, nY, zNet, nZ);
	DivideArea(xNet, nX, yNet, nY, zNet, nZ);

	Output();

	genT3D();
}
