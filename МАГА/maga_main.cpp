#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <set>
#include <locale>

using namespace std;

struct colour // ���� �����
{
	int red;
	int green;
	int blue;
};

struct point // �����
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
		return lhs.x == rhs.x
			&& lhs.y == rhs.y
			&& lhs.z == rhs.z;
	}

	friend bool operator!=(const point& lhs, const point& rhs)
	{
		return !(lhs == rhs);
	}
};

struct locateOfPoint
{
	int i, j, z;

	locateOfPoint(int i, int j, int z)
		: i(i),
		  j(j),
		  z(z)
	{
	}

	locateOfPoint(): i(0), j(0), z(0)
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

struct sigmStruct
{
	int terminalNode;
	vector<int> neighbors;
	vector<double> tNeighbors;
}*sigmT;
int *igT, *jgT;
double *ggT;

int LosLU(double* ggl, double* ggu, double* diag, int N, int* ig, int* jg, double* f, double* q);

vector<point> xyz_points;
vector<nvtr> KE;
vector<field> sreda;

double leftX, rightX, leftY, rightY, leftZ, rightZ;
double koordSourceX, koordSourceY, koordSourceZ;
double *xNet, *yNet, *zNet;
int nX, nY, nZ;
char*** matrixNode;
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

void inputConfig()
{
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
			cout << "������ �������� ����� �� ���. ���������� �������� ������ ���� ������ ��� ����� 1";
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

	cout << "������. �� ������� �����." << endl;
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

void OptimizationQuarterXforDuplication(int directionX, int directionY, int startX, int startY, int endX, int endY)
{
	int i, j, k, granX, granY, posI, posJ, posK;
	double height_1, height_2, width_1, width_2;
	for (j = startY * directionY; j < endY * directionY; j++)
	{
		posJ = j * directionY;
		for (i = startX * directionX; i < endX * directionX - 1; i++)
		{
			posI = i * directionX;
			//��������� ������� ��� ������� ����������
			if (directionX == 1 && directionY == 1 && (matrixNode[posI][posJ][0] == 'S' || matrixNode[posI][posJ][0] == 'R' || matrixNode[posI][posJ][0] == 'A' || matrixNode[posI][posJ][0] == 'V')
				|| (directionX == -1 && directionY == 1 && (matrixNode[posI][posJ][0] == 'S' || matrixNode[posI][posJ][0] == 'R' || matrixNode[posI][posJ][0] == 'D' || matrixNode[posI][posJ][0] == 'V'))
				|| (directionX == -1 && directionY == -1 && (matrixNode[posI][posJ][0] == 'W' || matrixNode[posI][posJ][0] == 'R' || matrixNode[posI][posJ][0] == 'D' || matrixNode[posI][posJ][0] == 'V'))
				|| (directionX == 1 && directionY == -1 && (matrixNode[posI][posJ][0] == 'W' || matrixNode[posI][posJ][0] == 'R' || matrixNode[posI][posJ][0] == 'A' || matrixNode[posI][posJ][0] == 'V')))
			{
				height_1 = height_2 = width_1 = width_2 = granX = granY = 0;
				for (k = j + 1; k <= endY * directionY; k++) //������ 1
				{
					posK = k * directionY;
					if (directionX == 1 && directionY == 1 && (matrixNode[posI][posK][0] == 'W' || matrixNode[posI][posK][0] == 'A' || matrixNode[posI][posK][0] == 'R' || matrixNode[posI][posK][0] == 'V')
						|| (directionX == -1 && directionY == 1 && (matrixNode[posI][posK][0] == 'W' || matrixNode[posI][posK][0] == 'D' || matrixNode[posI][posK][0] == 'R' || matrixNode[posI][posK][0] == 'V'))
						|| (directionX == -1 && directionY == -1 && (matrixNode[posI][posK][0] == 'S' || matrixNode[posI][posK][0] == 'D' || matrixNode[posI][posK][0] == 'R' || matrixNode[posI][posK][0] == 'V'))
						|| (directionX == 1 && directionY == -1 && (matrixNode[posI][posK][0] == 'S' || matrixNode[posI][posK][0] == 'A' || matrixNode[posI][posK][0] == 'R' || matrixNode[posI][posK][0] == 'V')))
					{
						height_1 = fabs(yNet[posK] - yNet[posJ]);
						break;
					}
				}
				for (k = i + 1; k < endX * directionX; k++) //������ 1
				{
					posK = k * directionX;
					if (directionY == 1 && (matrixNode[posK][posJ][0] == 'S' || matrixNode[posK][posJ][0] == 'R') ||
						directionY == -1 && (matrixNode[posK][posJ][0] == 'W' || matrixNode[posK][posJ][0] == 'R'))
					{
						width_1 = fabs(xNet[posK] - xNet[posI]);
						granX = posK;
						break;
					}
				}
				for (k = j + 1; k <= endY * directionY; k++) //������ 2
				{
					posK = k * directionY;
					if (directionY == 1 && (matrixNode[granX][posK][0] == 'W' || matrixNode[granX][posK][0] == 'R') ||
						directionY == -1 && (matrixNode[granX][posK][0] == 'S' || matrixNode[granX][posK][0] == 'R'))
					{
						height_2 = fabs(yNet[posK] - yNet[posJ]);
						granY = posK;
						break;
					}
				}
				if (height_2 == height_1)
				{
					for (k = granX * directionX + 1; k <= endX * directionX; k++)
					{
						posK = k * directionX;
						if (directionX == 1 && directionY == 1 && (matrixNode[posK][posJ][0] == 'D' || matrixNode[posK][posJ][0] == 'S' || matrixNode[posK][posJ][0] == 'R' || matrixNode[posK][posJ][0] == 'V') ||
							directionX == -1 && directionY == 1 && (matrixNode[posK][posJ][0] == 'A' || matrixNode[posK][posJ][0] == 'S' || matrixNode[posK][posJ][0] == 'R' || matrixNode[posK][posJ][0] == 'V') ||
							directionX == -1 && directionY == -1 && (matrixNode[posK][posJ][0] == 'A' || matrixNode[posK][posJ][0] == 'W' || matrixNode[posK][posJ][0] == 'R' || matrixNode[posK][posJ][0] == 'V') ||
							directionX == 1 && directionY == -1 && (matrixNode[posK][posJ][0] == 'D' || matrixNode[posK][posJ][0] == 'W' || matrixNode[posK][posJ][0] == 'R' || matrixNode[posK][posJ][0] == 'V'))
						{
							width_2 = width_1 + fabs(xNet[posK] - xNet[granX]);
							break;
						}
					}
					if (LikeASquare(width_1, height_1) > LikeASquare(width_2, height_2))
					{
						if (directionY == 1)
						{
							if (matrixNode[granX][posJ][0] == 'S')
								matrixNode[granX][posJ][0] = 'Y';
							else
								matrixNode[granX][posJ][0] = 'W';
							if (matrixNode[granX][granY][0] == 'W')
								matrixNode[granX][granY][0] = 'Y';
							else
								matrixNode[granX][granY][0] = 'S';
						}
						else
						{
							if (matrixNode[granX][posJ][0] == 'W')
								matrixNode[granX][posJ][0] = 'Y';
							else
								matrixNode[granX][posJ][0] = 'S';
							if (matrixNode[granX][granY][0] == 'S')
								matrixNode[granX][granY][0] = 'Y';
							else
								matrixNode[granX][granY][0] = 'W';
						}
					}
				}
			}
		}
	}
}

void normDividing()
{
	for (int i = 0; i < nX; i++)
	{
		matrixNode[i] = new char*[nY];
		for (int j = 0; j < nY; j++)
		{
			matrixNode[i][j] = new char[nZ];
			for (int k = 0; k < nZ; k++)
			{
				if (i == 0)
					matrixNode[i][j][k] = 'A';
				else if (i == nX - 1)
					matrixNode[i][j][k] = 'D';
				else if (j == 0)
					matrixNode[i][j][k] = 'S';
				else if (j == nY - 1)
					matrixNode[i][j][k] = 'W';
				else if (k == 0)
					matrixNode[i][j][k] = 'F';
				else if (k == nZ - 1)
					matrixNode[i][j][k] = 'B';
				else
					matrixNode[i][j][k] = 'R';
				if ((k == 0 || k == nZ - 1) && (i == 0 || j == 0 || i == nX - 1 || j == nY - 1) ||
					(i == 0 || i == nX - 1) && (k == 0 || j == 0 || k == nZ - 1 || j == nY - 1) ||
					(j == 0 || j == nY - 1) && (i == 0 || k == 0 || i == nX - 1 || k == nZ - 1))
					matrixNode[i][j][k] = 'V';
			}
		}
	}
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

void duplicatingOfXY()
{
	for (int i = 0; i < nX; i++)
	{
		matrixNode[i] = new char*[nY];
		for (int j = 0; j < nY; j++)
		{
			matrixNode[i][j] = new char[nZ];
			int k = 0;

			if (i == 0)
				matrixNode[i][j][k] = 'A';
			else if (i == nX - 1)
				matrixNode[i][j][k] = 'D';
			else if (j == 0)
				matrixNode[i][j][k] = 'S';
			else if (j == nY - 1)
				matrixNode[i][j][k] = 'W';
			else
				matrixNode[i][j][k] = 'R';
		}
	}
	matrixNode[0][0][0] = matrixNode[0][nY - 1][0] = matrixNode[nX - 1][0][0] = matrixNode[nX - 1][nY - 1][0] = 'V';

	double height_1, height_2, width_1, width_2;
	int locateSourceX, locateSourceY;
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

		//1-� ��������
		OptimizationQuarterXforDuplication(1, 1, locateSourceX, locateSourceY, FindLocate(xNet, nX, sreda[indSreda].x2), FindLocate(yNet, nY, sreda[indSreda].y2));
		//OptimizationQuarterY(1, 1, locateSourceX, locateSourceY, FindLocate(xNet, nX, sreda[indSreda].x2), FindLocate(yNet, nY, sreda[indSreda].y2));

		//2-� ��������
		OptimizationQuarterXforDuplication(-1, 1, locateSourceX, locateSourceY, FindLocate(xNet, nX, sreda[indSreda].x1), FindLocate(yNet, nY, sreda[indSreda].y2));
		//OptimizationQuarterY(-1, 1, locateSourceX, locateSourceY, FindLocate(xNet, nX, sreda[indSreda].x1), FindLocate(yNet, nY, sreda[indSreda].y2));

		//3-� ��������
		OptimizationQuarterXforDuplication(-1, -1, locateSourceX, locateSourceY, FindLocate(xNet, nX, sreda[indSreda].x1), FindLocate(yNet, nY, sreda[indSreda].y1));
		//OptimizationQuarterY(-1, -1, locateSourceX, locateSourceY, FindLocate(xNet, nX, sreda[indSreda].x1), FindLocate(yNet, nY, sreda[indSreda].y1));

		//4-� ��������
		OptimizationQuarterXforDuplication(1, -1, locateSourceX, locateSourceY, FindLocate(xNet, nX, sreda[indSreda].x2), FindLocate(yNet, nY, sreda[indSreda].y1));
		//OptimizationQuarterY(1, -1, locateSourceX, locateSourceY, FindLocate(xNet, nX, sreda[indSreda].x2), FindLocate(yNet, nY, sreda[indSreda].y1));
	}

	for (int k = 1; k < nZ; k++)
	{
		for (int i = 0; i < nX; i++)
			for (int j = 0; j < nY; j++)
			{
				matrixNode[i][j][k] = matrixNode[i][j][0];
			}
	}
	for (int i = 0; i < nX;i++)
	{
		for (int j = 0; j < nY;j++)
		{
			if (i == 0 || i == nX - 1 || j == 0 || j == nY - 1)
			{
				if(matrixNode[i][j][0]!='Y' && matrixNode[i][j][0] != 'V')
				{
					matrixNode[i][j][0] = matrixNode[i][j][nZ - 1] = 'V';
				}
			} else
			{
				if(matrixNode[i][j][0]=='R')
				{
					matrixNode[i][j][0] = 'F';
					matrixNode[i][j][nZ - 1] = 'B';
				}
			}
		}
	}
}

//W - ����� ��� �����
//A - ����� ��� �����
//S - ������� ��� �����
//D - ������ ��� �����
//F - ����� ��� �����
//B - ������ ��� �����

//V - �������
//R - ���������� ����, 6 ����� ���������
//Y - ����������� ������� ��� ����� �� �����
void DivideArea(double* xNet, int nX, double* yNet, int nY, double* zNet, int nZ)
{
	matrixNode = new char**[nX];
	for (int i = 0; i < nX; i++)
	{
		matrixNode[i] = new char*[nY];
		for (int j = 0; j < nY; j++)
		{
			matrixNode[i][j] = new char[nZ];

			/*for (int k = 0; k < nZ; k++)
			{
				if (i == 0)
					matrixNode[i][j][k] = 'A';
				else if (i == nX - 1)
					matrixNode[i][j][k] = 'D';
				else if (j == 0)
					matrixNode[i][j][k] = 'S';
				else if (j == nY - 1)
					matrixNode[i][j][k] = 'W';
				else if (k == 0)
					matrixNode[i][j][k] = 'F';
				else if (k == nZ - 1)
					matrixNode[i][j][k] = 'B';
				else
					matrixNode[i][j][k] = 'R';
				if ((k == 0 || k == nZ - 1) && (i == 0 || j == 0 || i == nX - 1 || j == nY - 1) ||
					(i == 0 || i == nX - 1) && (k == 0 || j == 0 || k == nZ - 1 || j == nY - 1) ||
					(j == 0 || j == nY - 1) && (i == 0 || k == 0 || i == nX - 1 || k == nZ - 1))
					matrixNode[i][j][k] = 'V';
			}*/
		}
	}
	duplicatingOfXY();

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
		cout << "������. �� ��������� ������ ����������. ����� ������� �� ����� ����� �������� �����������." << endl;
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
				// ���� ��� ���� �� �������� �������������
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
		//else
		//	for (k = igT[KE[ielem].uzel[i] - kolvoRegularNode]; k < igT[KE[ielem].uzel[i] - kolvoRegularNode + 1]; k++)
		//	{
		//		b[jgT[k]] += localB[i] * ggT[k];
		//	}
		for (j = 0; j < 8; j++)
		{
			if (KE[ielem].uzel[i] < kolvoRegularNode)
			{
				if (KE[ielem].uzel[j] < kolvoRegularNode)
				{
					AddToMatrix(KE[ielem].uzel[i], KE[ielem].uzel[j], localMatrix[i][j]);
				}
				/*else
				{
					for (k = igT[KE[ielem].uzel[j] - kolvoRegularNode]; k < igT[KE[ielem].uzel[j] - kolvoRegularNode + 1]; k++)
					{
						posJ = jgT[k];
						koefJ = ggT[k];
						AddToMatrix(KE[ielem].uzel[i], posJ, koefJ*localMatrix[i][j]);
					}
				}
			}
			else
			{
				if (KE[ielem].uzel[j] < kolvoRegularNode)
				{
					for (k = igT[KE[ielem].uzel[i] - kolvoRegularNode]; k < igT[KE[ielem].uzel[i] - kolvoRegularNode + 1]; k++)
					{
						posI = jgT[k];
						koefI = ggT[k];
						AddToMatrix(posI, KE[ielem].uzel[j], koefI*localMatrix[i][j]);
					}
				}
				else
				{
					for (k = igT[KE[ielem].uzel[i] - kolvoRegularNode]; k < igT[KE[ielem].uzel[i] - kolvoRegularNode + 1]; k++)
					{
						for (l = igT[KE[ielem].uzel[j] - kolvoRegularNode]; l < igT[KE[ielem].uzel[j] - kolvoRegularNode + 1]; l++)
						{
							AddToMatrix(jgT[k], jgT[l], ggT[k] * ggT[l] * localMatrix[i][j]);
						}
					}
				}*/
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

	cout << "������ findKE. �� ������ ��." << endl;
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

//	��������� ������� �� ������
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

	for (int k = 0; k<nZ; k++)
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
	cout << "������ � FindAreaNumber: �� ������� ����������." << endl;
	system("pause");
	exit(1);
	return -1;
}

void outputKEandXYZ()
{
	ofstream fileXY("xyz.txt");
	ofstream fileNvtr("nvtr.txt");
	vector<point> ncPoint;

	//��������� ������� ���������� � ������������ ������
	for (int k = 0; k < nZ;k++)
	for (int j = 0; j < nY; j++)
		for (int i = 0; i < nX; i++)
		{
			if (matrixNode[i][j][k] != 'Y')
				if ((k == 0 && matrixNode[i][j][k]=='F' || k == nZ - 1 && matrixNode[i][j][k]=='B'
					|| i == 0 && matrixNode[i][j][k]=='A' || i == nX - 1 && matrixNode[i][j][k]=='D'
					|| j == 0 && matrixNode[i][j][k] =='S' || j == nY - 1 && matrixNode[i][j][k]=='W') 
					|| matrixNode[i][j][k]=='V' || matrixNode[i][j][k] == 'R')
				{
					xyz_points.push_back(point(xNet[i], yNet[j], zNet[k]));
				}
				else
				{
					ncPoint.push_back(point(xNet[i], yNet[j], zNet[k]));
				}
		}

	//��������� � ����� ������� ������ ������������ �������
	for (vector<point>::iterator it = ncPoint.begin(); it < ncPoint.end(); it++)
	{
		xyz_points.push_back(*it);
	}

	//���������� �������� � T
	nColT = ncPoint.size();
	kolvoRegularNode = xyz_points.size() - nColT;

	//��������� ���� xy.txt
	fileXY << xyz_points.size() << " " << kolvoRegularNode << endl;
	for (int i = 0; i < xyz_points.size(); i++)
	{
		fileXY << xyz_points[i].x << " " << xyz_points[i].y << " " << xyz_points[i].z << endl;
	}

	//��������� ��������� ��
	nvtr tempNvtr;
	for (int iZ = 0; iZ < nZ - 1; iZ++)
	for (int iY = 0; iY < nY - 1; iY++)
		for (int iX = 0; iX < nX - 1; iX++)
		{
			if (matrixNode[iX][iY][iZ] == 'V' || matrixNode[iX][iY][iZ] == 'S' || matrixNode[iX][iY][iZ] == 'A'
				|| matrixNode[iX][iY][iZ] == 'R' || matrixNode[iX][iY][iZ] == 'F')
			{
				tempNvtr.uzel[0] = indexXYZ(xNet[iX], yNet[iY],zNet[iZ]);
				for (int k = iX + 1; k < nX; k++)
				{
					if (matrixNode[k][iY][iZ] == 'S' || matrixNode[k][iY][iZ] == 'R' || matrixNode[k][iY][iZ] == 'D' 
						|| matrixNode[k][iY][iZ] == 'V' || matrixNode[k][iY][iZ] == 'F')
					{
						tempNvtr.uzel[1] = indexXYZ(xNet[k], yNet[iY],zNet[iZ]);
						for (int t = iY + 1; t < nY; t++)
						{
							if (matrixNode[k][t][iZ] == 'W' || matrixNode[k][t][iZ] == 'R' || matrixNode[k][t][iZ] == 'D' 
								|| matrixNode[k][t][iZ] == 'V' || matrixNode[k][t][iZ] == 'F')
							{
								tempNvtr.uzel[3] = indexXYZ(xNet[k], yNet[t],zNet[iZ]);
								tempNvtr.uzel[2] = indexXYZ(xNet[iX], yNet[t], zNet[iZ]);

								for (int h = iZ; h < nZ;h++)
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

	//��������� ���� nvtr.txt
	fileNvtr << KE.size() << endl;
	for (int i = 0; i < KE.size(); i++)
	{
		for (int j = 0; j < 8; j++)
		{
			fileNvtr << KE[i].uzel[j] << " ";
		}
		fileNvtr <<KE[i].numberField<< endl;
	}
}

//void GenerateT()
//{
//	ofstream outT("Tmatrix.txt");
//	igT = new int[nColT + 1];
//
//	//��������� ��������������� ����-���������
//	int i, j, k, indNeib1, indNeib2;
//	sigmT = new sigmStruct[kolvoRegularNode];
//	for (i = kolvoRegularNode; i < xyz_points.size(); i++)
//	{
//		sigmT[i - kolvoRegularNode].terminalNode = i;
//		locateOfPoint term = FindLocate(xyz_points[i]);
//		if (matrixNode[term.i][term.j] == 'W' || matrixNode[term.i][term.j] == 'S')
//		{
//			for (indNeib1 = term.i - 1; indNeib1 >= 0 && matrixNode[indNeib1][term.j] != 'A' && matrixNode[indNeib1][term.j] != 'R'; indNeib1--);
//			sigmT[i - kolvoRegularNode].neighbors.push_back(NumberNode(xNet[indNeib1], yNet[term.j]));
//			for (indNeib2 = term.i + 1; indNeib2 < nX && matrixNode[indNeib2][term.j] != 'D' && matrixNode[indNeib2][term.j] != 'R'; indNeib2++);
//			sigmT[i - kolvoRegularNode].neighbors.push_back(NumberNode(xNet[indNeib2], yNet[term.j]));
//			sigmT[i - kolvoRegularNode].tNeighbors.push_back(fabs((xNet[indNeib2] - xNet[term.i]) / (xNet[indNeib2] - xNet[indNeib1])));
//			sigmT[i - kolvoRegularNode].tNeighbors.push_back(fabs((xNet[indNeib1] - xNet[term.i]) / (xNet[indNeib2] - xNet[indNeib1])));
//		}
//		else
//		{
//			for (indNeib1 = term.j - 1; indNeib1 >= 0 && matrixNode[term.i][indNeib1] != 'S' && matrixNode[term.i][indNeib1] != 'R'; indNeib1--);
//			sigmT[i - kolvoRegularNode].neighbors.push_back(NumberNode(xNet[term.i], yNet[indNeib1]));
//			for (indNeib2 = term.j + 1; indNeib2 < nY && matrixNode[term.i][indNeib2] != 'W' && matrixNode[term.i][indNeib2] != 'R'; indNeib2++);
//			sigmT[i - kolvoRegularNode].neighbors.push_back(NumberNode(xNet[term.i], yNet[indNeib2]));
//			sigmT[i - kolvoRegularNode].tNeighbors.push_back(fabs((yNet[indNeib2] - yNet[term.j]) / (yNet[indNeib2] - yNet[indNeib1])));
//			sigmT[i - kolvoRegularNode].tNeighbors.push_back(fabs((yNet[indNeib1] - yNet[term.j]) / (yNet[indNeib2] - yNet[indNeib1])));
//		}
//	}
//
//	for (i = 0; i < nColT; i++)
//	{
//		for (j = 0; j < 2; j++)
//		{
//			if (sigmT[i].neighbors[j] >= kolvoRegularNode)
//			{
//				sigmTChain(i, sigmT[i].neighbors[j] - kolvoRegularNode, sigmT[i].tNeighbors[j]);
//			}
//		}
//	}
//
//	for (i = 0; i < nColT; i++)
//		for (j = 1; j >= 0; j--)
//			if (sigmT[i].neighbors[j] >= kolvoRegularNode)
//			{
//				sigmT[i].neighbors.erase(sigmT[i].neighbors.begin() + j);
//				sigmT[i].tNeighbors.erase(sigmT[i].tNeighbors.begin() + j);
//			}
//
//	//��������� �� �����
//	double tmp_sort;
//	for (k = 0; k < nColT; k++)
//		for (i = 0; i < sigmT[k].neighbors.size(); ++i) // i - ����� �������� ����
//		{
//			int pos = i;
//			tmp_sort = sigmT[k].neighbors[i];
//			for (int j = i + 1; j < sigmT[k].neighbors.size(); ++j) // ���� ������ ����������� ��������
//			{
//				if (sigmT[k].neighbors[j] < tmp_sort)
//				{
//					pos = j;
//					tmp_sort = sigmT[k].neighbors[j];
//				}
//			}
//			sigmT[k].neighbors[pos] = sigmT[k].neighbors[i];
//			sigmT[k].neighbors[i] = tmp_sort; // ������ ������� ���������� � a[i]
//
//			tmp_sort = sigmT[k].tNeighbors[pos];
//			sigmT[k].tNeighbors[pos] = sigmT[k].tNeighbors[i];
//			sigmT[k].tNeighbors[i] = tmp_sort; // ������ ������� ���������� � a[i]
//
//		}
//
//	for (i = 0; i < nColT; i++)
//	{
//		outT << sigmT[i].terminalNode << " ";
//		for (j = 0; j < sigmT[i].neighbors.size(); j++)
//		{
//			outT << sigmT[i].neighbors[j] << " " << sigmT[i].tNeighbors[j] << " ";
//		}
//		outT << endl;
//	}
//	igT[0] = 0;
//	for (i = 0; i < nColT; i++)
//	{
//		igT[i + 1] = igT[i] + sigmT[i].neighbors.size();
//	}
//	ggT = new double[igT[nColT]];
//	jgT = new int[igT[nColT]];
//	for (i = 0, j = 0; i < nColT && j < igT[nColT]; i++)
//	{
//		for (k = 0; k < igT[i + 1] - igT[i]; k++, j++)
//		{
//			jgT[j] = sigmT[i].neighbors[k];
//			ggT[j] = sigmT[i].tNeighbors[k];
//		}
//	}
//
//	outT << endl << "T" << endl;
//	for (i = 0; i <= nColT; i++)
//	{
//		outT << igT[i] << " ";
//	}
//	outT << endl;
//	for (i = 0; i < igT[nColT]; i++)
//	{
//		outT << jgT[i] << " ";
//	}
//	outT << endl;
//	for (i = 0; i < igT[nColT]; i++)
//	{
//		outT << ggT[i] << " ";
//	}
//	outT << endl;
//
//}

int main(int argc, char* argv[])
{
	inputConfig();
	inputNet();
	DivideArea(xNet, nX, yNet, nY, zNet, nZ);
	//outputWithoutOptimisation();
	outputKEandXYZ();
	/*generatePortrait();
	GenerateMatrix();

	bool flagLU = false;
	if (flagLU)
	{
		int numberIteration = LosLU(ggl, ggu, di, xyz_points.size(), ig, jg, b, q);
		output << "numberIteration = " << numberIteration << endl;
	}
	else
	{
		runLOS();
	}

	calcPogreshnost(output);*/
}
