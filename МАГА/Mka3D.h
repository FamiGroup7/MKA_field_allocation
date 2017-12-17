#pragma once
#include "MkaUtils.h"
#include "Point.h"
#include "Edge.h"
#include "Qube.h"
//#include "LosLU.cpp"

using namespace std;

class Mka3D {
public:
	Mka3D();
	~Mka3D();

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

	struct locateOfPoint
	{
		int i, j, k;

		locateOfPoint(int i, int j, int k)
			: i(i),
			j(j),
			k(k)
		{
		}

		locateOfPoint() : i(0), j(0), k(0)
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

private:
	ofstream output;
	ofstream logger;
	const bool DEBUG = true;
	bool GRID_UNION = false;

	static const size_t BIT_SIZE = 13;
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

	const byte X_EDGE_DELETED = 0;
	const byte Y_EDGE_DELETED = 1;
	const byte Z_EDGE_DELETED = 2;

	int *igT, *jgT;
	double* ggT;
	vector<Point> xyz_points;
	vector<nvtr> KE;
	vector<field> sreda;
	multimap<Point, byte> termNodeOnEdge;
	double leftX, rightX, leftY, rightY, leftZ, rightZ;
	double koordSourceX, koordSourceY, koordSourceZ;
	double *xNet, *yNet, *zNet;
	int nX, nY, nZ;
	char*** matrixNode;
	bitset<Mka3D::BIT_SIZE>*** newNodes;
	int nColT, kolvoRegularNode;
	const int AXIS_SIZE = 6;
	enum Axis { LEFT, RIGHT, DOWN, UP, BACK, FORE };
	vector<sigmStruct3D> sigmNewT;

	int *ig, *jg;
	double *di, *b, *q, *ggl, *ggu;
	double localB[8], localMatrix[8][8];
	double G1[2][2] = { { 1,-1 },{ -1,1 } };
	double M1[2][2] = { { 1 / 3.,1 / 6. },{ 1 / 6.,1 / 3. } };
	double M2[4][4] = {
		{ 4,2,2,1 },
		{ 2,4,1,2 },
		{ 2,1,4,2 },
		{ 1,2,2,4 }
	};

	void inputConfig(); 
	void logError(char* message);
	void GenerateNetLikeTelma(set<double>& mas, ifstream & fileNet);
	int indexXYZ(Point goal);
	int indexXYZ(double x, double y, double z);
	double LikeASquare(double x, double y);
	double otn(double val1, double val2);
	double LikeACube(double x, double y, double z);
	int FindLocate(double * massiv, int razm, double x);
	int FindAreaNumber(int nodes[]);
	bool hasLeft(bitset<BIT_SIZE> node);
	bool hasRight(bitset<BIT_SIZE> node);
	bool hasUp(bitset<BIT_SIZE> node);
	bool hasDown(bitset<BIT_SIZE> node);
	bool hasBack(bitset<BIT_SIZE> node);
	bool hasFore(bitset<BIT_SIZE> node);
	bool hasXLines(bitset<BIT_SIZE> node);
	bool hasYLines(bitset<BIT_SIZE> node);
	bool hasZLines(bitset<BIT_SIZE> node);
	bool hasAll(bitset<BIT_SIZE> node);
	vector<Axis>* nodeInfo(bitset<BIT_SIZE> node);
	bool canOptimize(int directionX, int directionY, int directionZ, bitset<BIT_SIZE> node);
	int findNextX(int directionX, int directionY, int directionZ, int from_x, int to_x, int u_j, int u_t);
	int findNextY(int directionX, int directionY, int directionZ, int from_y, int to_y, int posI, int posT);
	int findNextZ(int directionX, int directionY, int directionZ, int from_z, int to_z, int posI, int posJ);
	void deletePlaneX(int xPlane, int x1, int x2, int y1, int y2, int z1, int z2);
	void deletePlaneY(int yPlane, int x1, int x2, int y1, int y2, int z1, int z2);
	Qube* getQube(int directionX, int directionY, int directionZ, int from_x, int from_y, int from_z, int to_x, int to_y, int to_z);
	void OptimizationQuarterX(int directionX, int directionY, int directionZ, int startX, int startY, int startZ, int endX, int endY, int endZ);
	void OptimizationQuarterY(int directionX, int directionY, int directionZ, int startX, int startY, int startZ, int endX, int endY, int endZ);
	void initNet(double * xNet, int nX, double * yNet, int nY, double * zNet, int nZ);
	void DivideArea(double * xNet, int nX, double * yNet, int nY, double * zNet, int nZ);
	void PrintLocalMatrix();
	void PrintPlotMatrix(bool flag_simmeric);
	void inputNet();
	void generatePortraitNesoglas();
	double analiticSolution(Point goal);
	double analiticSolution(double x, double y, double z);
	double Lambda(int ielem);
	double Gamma(int ielem);
	double Func(int ielem, int node);
	void CreateLocalMatrixs(int ielem);
	void AddToMatrix(int posI, int posJ, double el);
	void AddToA(int i, int j, double value);
	void AddToB(int i, double value);
	void Addition(int ielem);
	void doEdge1(ofstream & outEdge1File, const int intXorYorZ, const int kolvoRegularNode, double * varNet1, const int nVarNet1, double * varNet2, const int nVarNet2, const int unknownIndex);
	void Edge1_not_sim(bool up, bool down, bool left, bool right, bool fore, bool behind);
	int findKE(int ind_nodes[4]);
	void doEdge2(ofstream & outEdge2File, const int intXorYorZ, const int normalDirect, const int kolvoRegularNode, double * varNet1, const int nVarNet1, double * varNet2, const int nVarNet2, const int unknownIndex);
	void Edge2_not_sim(bool up, bool down, bool left, bool right, bool fore, bool behind);
	void doEdge3(ofstream & outEdge3File, int intXorYorZ, int normalDirect, int kolvoRegularNode, double * varNet1, int nVarNet1, double * varNet2, int nVarNet2, int unknownIndex);
	void Edge3_not_sim(bool up, bool down, bool left, bool right, bool fore, bool behind);
	void GenerateMatrix();
	void mult(double * res, double * v);
	double ScalarMult(double * v1, double * v2);
	void MultMatrixOnVector(double * in, double * out);
	void calcPogreshnost(ofstream & output);
	void runLOS();
	void outputWithoutOptimisation();
	locateOfPoint FindLocate(Point sample);
	void sigmTChain(int nTermNode, int startOfChain, double mnojT, set<int>& visitedNodes);
	set<byte> getValues(pair<map<Point, byte>::iterator, map<Point, byte>::iterator>& range);
	void genT3D();
	bool belongToInterval(double val, double l1, double l2);
	bool containsPointOnVergeOfKe(int termNode, int iKe);
	void constructXyzAndNvtr();
	void checkNet();

	//int LosLU(double* ggl, double* ggu, double* diag, int N, int* ig, int* jg, double* f, double* q);
};