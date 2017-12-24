#include "Mka3D.h"

int LosLU(double* ggl, double* ggu, double* diag, int N, int* ig, int* jg, double* f, double* q);

Mka3D::Mka3D(bool netOptimization, bool debugMod, bool optOnlyOnOneDirection, bool maxOptimization, bool optX, bool optY, bool optZ) :
	Mka3D("resources3D/", netOptimization, debugMod, optOnlyOnOneDirection, maxOptimization, optX, optY, optZ)
{
}

Mka3D::Mka3D(string filePrefix, bool netOptimization, bool debugMod, bool optOnlyOnOneDirection, bool maxOptimization,
	bool optX, bool optY, bool optZ)
{
	this->filePrefix = filePrefix;
	this->GRID_UNION = netOptimization;
	this->DEBUG = debugMod;
	this->X = optX;
	this->Y = optY;
	this->Z = optZ;
	this->optOnlyOnOneDirection = optOnlyOnOneDirection;
	this->MAX_OPTIMIZATION = maxOptimization;
	output.open(filePrefix + "solution.txt");
	logger.open(filePrefix + "mka3Dlog.txt");
	profiler.open(filePrefix + "profiling.txt");
	inputConfig();
}

void Mka3D::startFullProcess() {
	auto startProgramTime = std::chrono::system_clock::now();

	buildNet("sreda.txt", "sourceLocate.txt");

	netOptimization();

	prepareNetForSolve();

	build_xyz_nvtr_portratin_Abqx();
	generateGlobalMatrix();

	//LosLU(ggl, ggu, di, countRegularNodes, ig, jg, b, q);
	runLOS();
	calcPogreshnost(output);


	profiler << setw(40) << std::left << "Total program duration " <<
		MkaUtils::formattingTime(std::chrono::system_clock::time_point(std::chrono::system_clock::now() - startProgramTime)) << endl;

	profiler << endl;
	profiler << "KE number = " << KE.size() << endl;
	profiler << "dim(A) = " << countRegularNodes << endl;
	profiler << "Terminal nodes = " << nColT << endl;
}

void Mka3D::buildNet(string sredaFile, string sourceFile) {
	inputNet(filePrefix + sredaFile, filePrefix + sourceFile);
	initNet(xNet, nX, yNet, nY, zNet, nZ);
}

void Mka3D::build_xyz_nvtr_portratin_Abqx() {
	auto start = std::chrono::system_clock::now();
	constructXyzAndNvtr();
	profiler << setw(40) << std::left << "Construct xyz and nvtr " <<
		MkaUtils::formattingTime(std::chrono::system_clock::time_point(std::chrono::system_clock::now() - start)) << endl;

	start = std::chrono::system_clock::now();
	genT3D();
	profiler << setw(40) << std::left << "Generate T matrix " <<
		MkaUtils::formattingTime(std::chrono::system_clock::time_point(std::chrono::system_clock::now() - start)) << endl;

	start = std::chrono::system_clock::now();
	generatePortraitNesoglas();
	profiler << setw(40) << std::left << "Generate portrait " <<
		MkaUtils::formattingTime(std::chrono::system_clock::time_point(std::chrono::system_clock::now() - start)) << endl;
}

Mka3D::~Mka3D()
{
}

void Mka3D::inputConfig()
{
}

void Mka3D::logError(char* message) {
	cout << message << endl;
	throw exception(message);
}

void Mka3D::GenerateNetLikeTelma(set<double>& mas, ifstream& fileNet)
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
			cout << "ќшибка описани€ сетки по оси.  оэфициент разр€дки должен быть больше или равен 1";
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

int Mka3D::indexXYZ(Point goal)
{
	for (int i = 0; i < xyz_points.size(); i++)
	{
		if (goal == xyz_points[i])
			return i;
	}

	//cout << "ќшибка. Ќе найдена точка (" << goal.x << "," << goal.y << "," << goal.z << ")" << endl;
	//throw new exception();
	//system("pause");
	//exit(1);
	return -1;
}

int Mka3D::indexXYZ(double x, double y, double z)
{
	Point goal(x, y, z);
	return indexXYZ(goal);
}

double Mka3D::LikeASquare(double x, double y)
{
	if (x > y)
		return (1 - y / x);
	else
		return (1 - x / y);
}

double Mka3D::otn(double val1, double val2) {
	if (val1 > val2)return val2 / val1;
	return val1 / val2;
}

double Mka3D::LikeACube(double x, double y, double z) {
	//double xy = x*y;
	//double xz = x*z;
	//double yz = y*z;
	//return (otn(x, y) + otn(x, z) + otn(y, z)) / 3.0;
	double min = MkaUtils::compare(x, y) == -1 ? x : y;
	return MkaUtils::compare(min, z) == -1 ? min : z;
}

//TODO search in sorted array
int Mka3D::FindLocate(double* massiv, int razm, double x)
{
	int i;
	for (i = 0; i < razm; i++)
	{
		if (massiv[i] == x)
			return i;
	}
	return -1;
}

int Mka3D::FindAreaNumber(int nodes[])
{
	int i;
	for (i = 0; i < sreda.size(); i++)
	{
		if (xyz_points[nodes[0]].x >= sreda[i].x1 && xyz_points[nodes[1]].x <= sreda[i].x2
			&& xyz_points[nodes[0]].y >= sreda[i].y1 && xyz_points[nodes[2]].y <= sreda[i].y2
			&& xyz_points[nodes[0]].z >= sreda[i].z1 && xyz_points[nodes[4]].z <= sreda[i].z2)
			return i;
	}
	logError("ќшибка в FindAreaNumber: не найдена подобласть.");
	return -1;
}

bool Mka3D::hasLeft(bitset<BIT_SIZE> node) {
	return node.test(LEFT_BACK) || node.test(LEFT_FORE)
		|| node.test(LEFT_DOWN) || node.test(LEFT_UP);
}

bool Mka3D::hasRight(bitset<BIT_SIZE> node) {
	return node.test(RIGHT_BACK) || node.test(RIGHT_FORE)
		|| node.test(RIGHT_DOWN) || node.test(RIGHT_UP);
}

bool Mka3D::hasUp(bitset<BIT_SIZE> node) {
	return node.test(LEFT_UP) || node.test(RIGHT_UP)
		|| node.test(BACK_UP) || node.test(FORE_UP);
}

bool Mka3D::hasDown(bitset<BIT_SIZE> node) {
	return node.test(LEFT_DOWN) || node.test(RIGHT_DOWN)
		|| node.test(BACK_DOWN) || node.test(FORE_DOWN);
}

bool Mka3D::hasBack(bitset<BIT_SIZE> node) {
	return node.test(LEFT_BACK) || node.test(RIGHT_BACK)
		|| node.test(BACK_UP) || node.test(BACK_DOWN);
}

bool Mka3D::hasFore(bitset<BIT_SIZE> node) {
	return node.test(LEFT_FORE) || node.test(RIGHT_FORE)
		|| node.test(FORE_UP) || node.test(FORE_DOWN);
}

bool Mka3D::hasXLines(bitset<BIT_SIZE> node) {
	return hasLeft(node) && hasRight(node);
}

bool Mka3D::hasYLines(bitset<BIT_SIZE> node) {
	return hasBack(node) && hasFore(node);
}

bool Mka3D::hasZLines(bitset<BIT_SIZE> node) {
	return hasUp(node) && hasDown(node);
}

bool Mka3D::hasAll(bitset<BIT_SIZE> node) {
	return hasXLines(node) && hasYLines(node) && hasZLines(node);
}

vector<Mka3D::Axis>* Mka3D::nodeInfo(bitset<BIT_SIZE> node) {
	int n = 0;
	vector<Axis>* result = new vector<Axis>;
	if (hasLeft(node)) result->push_back(LEFT);
	if (hasRight(node)) result->push_back(RIGHT);
	if (hasDown(node)) result->push_back(DOWN);
	if (hasUp(node)) result->push_back(UP);
	if (hasBack(node)) result->push_back(BACK);
	if (hasFore(node)) result->push_back(FORE);
	return result;
}

bool Mka3D::canOptimize(int directionX, int directionY, int directionZ, bitset<BIT_SIZE> node) {
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

int Mka3D::findNextX(int directionX, int directionY, int directionZ,
	int from_x, int to_x, int u_j, int u_t) {
	int posK;
	for (int k = from_x*directionX + 1; k <= to_x * directionX; k++)
	{
		posK = k * directionX;
		if (canOptimize(-1 * directionX, directionY, directionZ, newNodes[posK][u_j][u_t]))
		{
			return posK;
		}
	}
	return -1;
}

int Mka3D::findNextY(int directionX, int directionY, int directionZ,
	int from_y, int to_y, int posI, int posT) {
	int posK;
	for (int k = from_y*directionY + 1; k <= to_y * directionY; k++)
	{
		posK = k * directionY;
		if (canOptimize(directionX, -1 * directionY, directionZ, newNodes[posI][posK][posT]))
		{
			return posK;
		}
	}
	return -1;
}

int Mka3D::findNextZ(int directionX, int directionY, int directionZ,
	int from_z, int to_z, int posI, int posJ) {
	int posK;
	for (int k = from_z*directionZ + 1; k <= to_z * directionZ; k++)
	{
		posK = k * directionZ;
		if (canOptimize(directionX, directionY, -1 * directionZ, newNodes[posI][posJ][posK]))
		{
			return posK;
		}
	}
	return -1;
}

bool Mka3D::testThatPointHasOptimizedOnAnotherDirection(int iX, int iY, int iZ, byte currentDirection) {
	pair<map<Point, byte>::iterator, map<Point, byte>::iterator> range
		= termNodeOnEdge.equal_range(Point(xNet[iX],yNet[iY],zNet[iZ]));
	set<byte> edges = getValues(range);
	if (edges.find(currentDirection) != edges.end()) {
		return edges.size() > 1;
	}
	else {
		return edges.size() > 0;
	}
}

bool Mka3D::deletePlaneX(int xPlane, int x1, int x2, int y1, int y2, int z1, int z2) {
	if(optOnlyOnOneDirection)
	if (
		testThatPointHasOptimizedOnAnotherDirection(xPlane, y1, z1, X_EDGE_DELETED) ||
		testThatPointHasOptimizedOnAnotherDirection(xPlane, y1, z2, X_EDGE_DELETED) ||
		testThatPointHasOptimizedOnAnotherDirection(xPlane, y2, z1, X_EDGE_DELETED) ||
		testThatPointHasOptimizedOnAnotherDirection(xPlane, y2, z2, X_EDGE_DELETED)
		)
	{
		logger << "Cannot delete plane x = " << xNet[xPlane] << ", y1 = " << yNet[y1] << ", y2 = " << yNet[y2]
			<< ", z1 = " << zNet[z1] << ", z2 = " << zNet[z2] << endl;
		return false;
	}

	logger << "Deleted plane x = " << xNet[xPlane] << ", y1 = " << yNet[y1] << ", y2 = " << yNet[y2]
		<< ", z1 = " << zNet[z1] << ", z2 = " << zNet[z2] << endl;

	newNodes[xPlane][y1][z1].reset(IS_REGULAR);
	newNodes[xPlane][y1][z2].reset(IS_REGULAR);
	newNodes[xPlane][y2][z1].reset(IS_REGULAR);
	newNodes[xPlane][y2][z2].reset(IS_REGULAR);

	newNodes[xPlane][y1][z1].reset(BACK_UP);
	newNodes[xPlane][y1][z2].reset(BACK_DOWN);
	newNodes[xPlane][y2][z1].reset(FORE_UP);
	newNodes[xPlane][y2][z2].reset(FORE_DOWN);

	//нижн€€
	if (!newNodes[xPlane][y1][z1].test(BACK_DOWN))
		newNodes[xPlane][y1][z1].reset(LEFT_BACK).reset(RIGHT_BACK);
	if (!newNodes[xPlane][y2][z1].test(FORE_DOWN))
		newNodes[xPlane][y2][z1].reset(LEFT_FORE).reset(RIGHT_FORE);
	//верхн€€
	if (!newNodes[xPlane][y1][z2].test(BACK_UP))
		newNodes[xPlane][y1][z2].reset(LEFT_BACK).reset(RIGHT_BACK);
	if (!newNodes[xPlane][y2][z2].test(FORE_UP))
		newNodes[xPlane][y2][z2].reset(LEFT_FORE).reset(RIGHT_FORE);
	//ближн€€
	if (!newNodes[xPlane][y1][z1].test(FORE_UP))
		newNodes[xPlane][y1][z1].reset(LEFT_UP).reset(RIGHT_UP);
	if (!newNodes[xPlane][y1][z2].test(FORE_DOWN))
		newNodes[xPlane][y1][z2].reset(LEFT_DOWN).reset(RIGHT_DOWN);
	//дальн€€
	if (!newNodes[xPlane][y2][z1].test(BACK_UP))
		newNodes[xPlane][y2][z1].reset(LEFT_UP).reset(RIGHT_UP);
	if (!newNodes[xPlane][y2][z2].test(BACK_DOWN))
		newNodes[xPlane][y2][z2].reset(LEFT_DOWN).reset(RIGHT_DOWN);

	termNodeOnEdge.insert(pair<Point, byte>(Point(xNet[xPlane], yNet[y1], zNet[z1]), X_EDGE_DELETED));
	termNodeOnEdge.insert(pair<Point, byte>(Point(xNet[xPlane], yNet[y1], zNet[z2]), X_EDGE_DELETED));
	termNodeOnEdge.insert(pair<Point, byte>(Point(xNet[xPlane], yNet[y2], zNet[z1]), X_EDGE_DELETED));
	termNodeOnEdge.insert(pair<Point, byte>(Point(xNet[xPlane], yNet[y2], zNet[z2]), X_EDGE_DELETED));
	return true;
}

bool Mka3D::deletePlaneY(int yPlane, int x1, int x2, int y1, int y2, int z1, int z2) {
	if (optOnlyOnOneDirection)
	if
		(
			testThatPointHasOptimizedOnAnotherDirection(x1, yPlane, z1, Y_EDGE_DELETED) ||
			testThatPointHasOptimizedOnAnotherDirection(x1, yPlane, z2, Y_EDGE_DELETED) ||
			testThatPointHasOptimizedOnAnotherDirection(x2, yPlane, z1, Y_EDGE_DELETED) ||
			testThatPointHasOptimizedOnAnotherDirection(x2, yPlane, z2, Y_EDGE_DELETED)
			)
	{
		logger << "Cannot delete plane y = " << yNet[yPlane] << ", x1 = " << xNet[x1] << ", x2 = " << xNet[x2]
			<< ", z1 = " << zNet[z1] << ", z2 = " << zNet[z2] << endl;
		return false;
	}

	logger << "Deleted plane y = " << yNet[yPlane] << ", x1 = " << xNet[x1] << ", x2 = " << xNet[x2]
		<< ", z1 = " << zNet[z1] << ", z2 = " << zNet[z2] << endl;

	newNodes[x1][yPlane][z1].reset(IS_REGULAR);
	newNodes[x1][yPlane][z2].reset(IS_REGULAR);
	newNodes[x2][yPlane][z1].reset(IS_REGULAR);
	newNodes[x2][yPlane][z2].reset(IS_REGULAR);

	newNodes[x1][yPlane][z1].reset(RIGHT_UP);
	newNodes[x1][yPlane][z2].reset(RIGHT_DOWN);
	newNodes[x2][yPlane][z1].reset(LEFT_UP);
	newNodes[x2][yPlane][z2].reset(LEFT_DOWN);

	//нижн€€
	if (!newNodes[x1][yPlane][z1].test(RIGHT_DOWN))
		newNodes[x1][yPlane][z1].reset(RIGHT_FORE).reset(RIGHT_BACK);
	if (!newNodes[x2][yPlane][z1].test(LEFT_DOWN))
		newNodes[x2][yPlane][z1].reset(LEFT_FORE).reset(LEFT_BACK);
	//верхн€€
	if (!newNodes[x1][yPlane][z2].test(RIGHT_UP))
		newNodes[x1][yPlane][z2].reset(RIGHT_FORE).reset(RIGHT_BACK);
	if (!newNodes[x2][yPlane][z2].test(LEFT_UP))
		newNodes[x2][yPlane][z2].reset(LEFT_FORE).reset(LEFT_BACK);
	//лева€
	if (!newNodes[x1][yPlane][z1].test(LEFT_UP))
		newNodes[x1][yPlane][z1].reset(BACK_UP).reset(FORE_UP);
	if (!newNodes[x1][yPlane][z2].test(LEFT_DOWN))
		newNodes[x1][yPlane][z2].reset(BACK_DOWN).reset(FORE_DOWN);
	//права€
	if (!newNodes[x2][yPlane][z1].test(RIGHT_UP))
		newNodes[x2][yPlane][z1].reset(BACK_UP).reset(FORE_UP);
	if (!newNodes[x2][yPlane][z2].test(RIGHT_DOWN))
		newNodes[x2][yPlane][z2].reset(BACK_DOWN).reset(FORE_DOWN);

	termNodeOnEdge.insert(pair<Point, byte>(Point(xNet[x1], yNet[yPlane], zNet[z1]), Y_EDGE_DELETED));
	termNodeOnEdge.insert(pair<Point, byte>(Point(xNet[x1], yNet[yPlane], zNet[z2]), Y_EDGE_DELETED));
	termNodeOnEdge.insert(pair<Point, byte>(Point(xNet[x2], yNet[yPlane], zNet[z1]), Y_EDGE_DELETED));
	termNodeOnEdge.insert(pair<Point, byte>(Point(xNet[x2], yNet[yPlane], zNet[z2]), Y_EDGE_DELETED));
	return true;
}

bool Mka3D::deletePlaneZ(int zPlane, int x1, int x2, int y1, int y2, int z1, int z2) {
	if (optOnlyOnOneDirection)
	if
		(
			testThatPointHasOptimizedOnAnotherDirection(x1, y1, zPlane, Z_EDGE_DELETED) ||
			testThatPointHasOptimizedOnAnotherDirection(x1, y2, zPlane, Z_EDGE_DELETED) ||
			testThatPointHasOptimizedOnAnotherDirection(x2, y1, zPlane, Z_EDGE_DELETED) ||
			testThatPointHasOptimizedOnAnotherDirection(x2, y2, zPlane, Z_EDGE_DELETED)
			)
	{
		logger << "Cannot delete plane z = " << zNet[zPlane] << ", x1 = " << xNet[x1] << ", x2 = " << xNet[x2]
			<< ", y1 = " << yNet[y1] << ", y2 = " << yNet[y2] << endl;
		return false;
	}

	logger << "Deleted plane z = " << zNet[zPlane] << ", x1 = " << xNet[x1] << ", x2 = " << xNet[x2]
		<< ", y1 = " << yNet[y1] << ", y2 = " << yNet[y2] << endl;

	newNodes[x1][y1][zPlane].reset(IS_REGULAR);
	newNodes[x1][y2][zPlane].reset(IS_REGULAR);
	newNodes[x2][y1][zPlane].reset(IS_REGULAR);
	newNodes[x2][y2][zPlane].reset(IS_REGULAR);

	newNodes[x1][y1][zPlane].reset(RIGHT_BACK);
	newNodes[x1][y2][zPlane].reset(RIGHT_FORE);
	newNodes[x2][y1][zPlane].reset(LEFT_BACK);
	newNodes[x2][y2][zPlane].reset(LEFT_FORE);

	//передн€€
	if (!newNodes[x1][y1][zPlane].test(RIGHT_FORE))
		newNodes[x1][y1][zPlane].reset(RIGHT_DOWN).reset(RIGHT_UP);
	if (!newNodes[x2][y1][zPlane].test(LEFT_FORE))
		newNodes[x2][y1][zPlane].reset(LEFT_DOWN).reset(LEFT_UP);
	//задн€€
	if (!newNodes[x1][y2][zPlane].test(RIGHT_BACK))
		newNodes[x1][y2][zPlane].reset(RIGHT_DOWN).reset(RIGHT_UP);
	if (!newNodes[x2][y2][zPlane].test(LEFT_BACK))
		newNodes[x2][y2][zPlane].reset(LEFT_DOWN).reset(LEFT_UP);
	//лева€
	if (!newNodes[x1][y1][zPlane].test(LEFT_BACK))
		newNodes[x1][y1][zPlane].reset(BACK_UP).reset(BACK_DOWN);
	if (!newNodes[x1][y2][zPlane].test(LEFT_FORE))
		newNodes[x1][y2][zPlane].reset(FORE_UP).reset(FORE_DOWN);
	//права€
	if (!newNodes[x2][y1][zPlane].test(RIGHT_BACK))
		newNodes[x2][y1][zPlane].reset(BACK_UP).reset(BACK_DOWN);
	if (!newNodes[x2][y2][zPlane].test(RIGHT_FORE))
		newNodes[x2][y2][zPlane].reset(FORE_UP).reset(FORE_DOWN);

	termNodeOnEdge.insert(pair<Point, byte>(Point(xNet[x1], yNet[y1], zNet[zPlane]), Z_EDGE_DELETED));
	termNodeOnEdge.insert(pair<Point, byte>(Point(xNet[x1], yNet[y2], zNet[zPlane]), Z_EDGE_DELETED));
	termNodeOnEdge.insert(pair<Point, byte>(Point(xNet[x2], yNet[y1], zNet[zPlane]), Z_EDGE_DELETED));
	termNodeOnEdge.insert(pair<Point, byte>(Point(xNet[x2], yNet[y2], zNet[zPlane]), Z_EDGE_DELETED));
	return true;
}

Qube* Mka3D::getQube(int directionX, int directionY, int directionZ,
	int from_x, int from_y, int from_z, int to_x, int to_y, int to_z)
{
	int i_nextX = findNextX(directionX, directionY, directionZ, from_x, to_x, from_y, from_z);
	int i_nextY = findNextY(directionX, directionY, directionZ, from_y, to_y, from_x, from_z);
	int i_nextZ = findNextZ(directionX, directionY, directionZ, from_z, to_z, from_x, from_y);

	if (i_nextX == -1 || i_nextY == -1 || i_nextZ == -1) {
		return NULL;
	}

	if (DEBUG) {
		int X_y1_z2 = findNextX(directionX, directionY, -1 * directionZ, from_x, to_x, from_y, i_nextZ);
		int X_y2_z1 = findNextX(directionX, -1 * directionY, directionZ, from_x, to_x, i_nextY, from_z);
		int X_y2_z2 = findNextX(directionX, -1 * directionY, -1 * directionZ, from_x, to_x, i_nextY, i_nextZ);
		if (X_y1_z2 == -1 || X_y1_z2 != X_y2_z1 || X_y2_z1 != X_y2_z2 || X_y2_z2 != i_nextX) {
			logError("Error in qube check X");
		}

		int Y_x1_z2 = findNextY(directionX, directionY, -1 * directionZ, from_y, to_y, from_x, i_nextZ);
		int Y_x2_z1 = findNextY(-1 * directionX, directionY, directionZ, from_y, to_y, i_nextX, from_z);
		int Y_x2_z2 = findNextY(-1 * directionX, directionY, -1 * directionZ, from_y, to_y, i_nextX, i_nextZ);
		if (Y_x1_z2 == -1 || Y_x1_z2 != Y_x2_z1 || Y_x2_z1 != Y_x2_z2 || Y_x2_z2 != i_nextY) {
			logError("Error in qube check Y");
		}

		int Z_x1_y2 = findNextZ(directionX, -1 * directionY, directionZ, from_z, to_z, from_x, i_nextY);
		int Z_x2_y1 = findNextZ(-1 * directionX, directionY, directionZ, from_z, to_z, i_nextX, from_y);
		int Z_x2_y2 = findNextZ(-1 * directionX, -1 * directionY, directionZ, from_z, to_z, i_nextX, i_nextY);
		if (Z_x1_y2 == -1 || Z_x1_y2 != Z_x2_y1 || Z_x2_y1 != Z_x2_y2 || Z_x2_y2 != i_nextZ) {
			logError("Error in qube check Z");
		}
	}
	return new Qube(from_x, from_y, from_z, i_nextX, i_nextY, i_nextZ);
}

void Mka3D::OptimizationQuarterX(int directionX, int directionY, int directionZ,
	int startX, int startY, int startZ, int endX, int endY, int endZ)
{
	int i, j, t, u_i, u_j, u_t;
	for (t = startZ * directionZ; t < endZ * directionZ; t++)
	{
		u_t = t * directionZ;
		for (j = startY * directionY; j < endY * directionY; j++)
		{
			u_j = j * directionY;
			for (i = startX * directionX; i < endX * directionX; i++)
			{
				u_i = i * directionX;
				if (canOptimize(directionX, directionY, directionZ, newNodes[u_i][u_j][u_t]))
				{
					Qube*current = getQube(directionX, directionY, directionZ, u_i, u_j, u_t, endX, endY, endZ);
					if (current == NULL) {
						continue;
					}
					if (!canOptimize(directionX, directionY, directionZ, newNodes[current->i_nextX][u_j][u_t])) {
						continue;
					}
					Qube*nextQube = getQube(directionX, directionY, directionZ, current->i_nextX, u_j, u_t, endX, endY, endZ);
					if (nextQube == NULL) {
						delete current;
						continue;
					}
					if (current->i_nextY == nextQube->i_nextY && current->i_nextZ == nextQube->i_nextZ) {
						double deep = current->getDepth(yNet);
						double height = current->getHeight(zNet);
						double width_1 = current->getWidth(xNet);
						double width_union = width_1 + nextQube->getWidth(xNet);

						if (LikeACube(width_1, deep, height) < LikeACube(width_union, deep, height)) {
							bool result;
							if (directionY == 1)
								if (directionZ == 1)
									result = deletePlaneX(current->i_nextX, u_i, nextQube->i_nextX, u_j, current->i_nextY, u_t, current->i_nextZ);
								else
									result = deletePlaneX(current->i_nextX, u_i, nextQube->i_nextX, u_j, current->i_nextY, current->i_nextZ, u_t);
							else
								if (directionZ == 1)
									result = deletePlaneX(current->i_nextX, u_i, nextQube->i_nextX, current->i_nextY, u_j, u_t, current->i_nextZ);
								else
									result = deletePlaneX(current->i_nextX, u_i, nextQube->i_nextX, current->i_nextY, u_j, current->i_nextZ, u_t);
							if (MAX_OPTIMIZATION && result) {
								i--;
							}
						}
					}
					delete current;
					delete nextQube;
				}
			}
		}
	}
}

void Mka3D::OptimizationQuarterY(int directionX, int directionY, int directionZ,
	int startX, int startY, int startZ, int endX, int endY, int endZ)
{
	int i, j, t, u_i, u_j, u_t;
	for (t = startZ * directionZ; t < endZ * directionZ; t++)
	{
		u_t = t * directionZ;
		for (i = startX * directionX; i < endX * directionX; i++)
		{
			u_i = i * directionX;
			for (j = startY * directionY; j < endY * directionY; j++)
			{
				u_j = j * directionY;
				if (canOptimize(directionX, directionY, directionZ, newNodes[u_i][u_j][u_t]))
				{
					Qube*current = getQube(directionX, directionY, directionZ, u_i, u_j, u_t, endX, endY, endZ);
					if (current == NULL) {
						continue;
					}
					if (!canOptimize(directionX, directionY, directionZ, newNodes[u_i][current->i_nextY][u_t])) {
						continue;
					}
					Qube*nextQube = getQube(directionX, directionY, directionZ, u_i, current->i_nextY, u_t, endX, endY, endZ);
					if (nextQube == NULL) {
						delete current;
						continue;
					}
					if (current->i_nextX == nextQube->i_nextX && current->i_nextZ == nextQube->i_nextZ) {
						double width = current->getWidth(xNet);
						double height = current->getHeight(zNet);
						double deep_1 = current->getDepth(yNet);
						double deep_union = deep_1 + nextQube->getDepth(yNet);

						if (LikeACube(width, deep_1, height) < LikeACube(width, deep_union, height)) {
							bool result;
							if (directionX == 1)
								if (directionZ == 1)
									result = deletePlaneY(current->i_nextY, u_i, current->i_nextX, u_j, nextQube->i_nextY, u_t, current->i_nextZ);
								else
									result = deletePlaneY(current->i_nextY, u_i, current->i_nextX, u_j, nextQube->i_nextY, current->i_nextZ, u_t);
							else
								if (directionZ == 1)
									result = deletePlaneY(current->i_nextY, current->i_nextX, u_i, u_j, nextQube->i_nextY, u_t, current->i_nextZ);
								else
									result = deletePlaneY(current->i_nextY, current->i_nextX, u_i, u_j, nextQube->i_nextY, current->i_nextZ, u_t);
							if (MAX_OPTIMIZATION && result) {
								j--;
							}
						}
					}
					delete current;
					delete nextQube;
				}
			}
		}
	}
}

void Mka3D::OptimizationQuarterZ(int directionX, int directionY, int directionZ,
	int startX, int startY, int startZ, int endX, int endY, int endZ)
{
	int i, j, t, u_i, u_j, u_t;
	for (i = startX * directionX; i < endX * directionX; i++)
	{
		u_i = i * directionX;
		for (j = startY * directionY; j < endY * directionY; j++)
		{
			u_j = j * directionY;
			for (t = startZ * directionZ; t < endZ * directionZ; t++)
			{
				u_t = t * directionZ;
				if (canOptimize(directionX, directionY, directionZ, newNodes[u_i][u_j][u_t]))
				{
					Qube*current = getQube(directionX, directionY, directionZ, u_i, u_j, u_t, endX, endY, endZ);
					if (current == NULL) {
						continue;
					}
					if (!canOptimize(directionX, directionY, directionZ, newNodes[u_i][u_j][current->i_nextZ])) {
						continue;
					}
					Qube*nextQube = getQube(directionX, directionY, directionZ, u_i, u_j, current->i_nextZ, endX, endY, endZ);
					if (nextQube == NULL) {
						delete current;
						continue;
					}
					if (current->i_nextX == nextQube->i_nextX && current->i_nextY == nextQube->i_nextY) {
						double width = current->getWidth(xNet);
						double deep = current->getDepth(yNet);
						double height_1 = current->getHeight(zNet);
						double height_union = height_1 + nextQube->getHeight(yNet);

						if (LikeACube(width, deep, height_1) < LikeACube(width, deep, height_union)) {
							bool result;
							if (directionX == 1)
								if (directionY == 1)
									result = deletePlaneZ(current->i_nextZ, u_i, current->i_nextX, u_j, current->i_nextY, u_t, nextQube->i_nextZ);
								else
									result = deletePlaneZ(current->i_nextZ, u_i, current->i_nextX, current->i_nextY, u_j, nextQube->i_nextZ, u_t);
							else
								if (directionY == 1)
									result = deletePlaneZ(current->i_nextZ, current->i_nextX, u_i, u_j, current->i_nextY, u_t, nextQube->i_nextZ);
								else
									result = deletePlaneZ(current->i_nextZ, current->i_nextX, u_i, current->i_nextY, u_j, nextQube->i_nextZ, u_t);
							if (MAX_OPTIMIZATION && result) {
								t--;
							}
						}
					}
					delete current;
					delete nextQube;
				}
			}
		}
	}
}

void Mka3D::initNet(double* xNet, int nX, double* yNet, int nY, double* zNet, int nZ) {
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

void Mka3D::netOptimization()
{
	double height_1, height_2, width_1, width_2;
	int locateSourceX, locateSourceY, locateSourceZ;
	if (!GRID_UNION) {
		logger << "Grid optimization disabled" << endl;
		return;
	}

	auto start = std::chrono::system_clock::now();
	for (int indSreda = 0; indSreda < sreda.size(); indSreda++)
	{
		int locateX1 = FindLocate(xNet, nX, sreda[indSreda].x1);
		int locateX2 = FindLocate(xNet, nX, sreda[indSreda].x2);
		int locateY1 = FindLocate(yNet, nY, sreda[indSreda].y1);
		int locateY2 = FindLocate(yNet, nY, sreda[indSreda].y2);
		int locateZ1 = FindLocate(zNet, nZ, sreda[indSreda].z1);
		int locateZ2 = FindLocate(zNet, nZ, sreda[indSreda].z2);

		if (koordSourceX <= sreda[indSreda].x1)
			locateSourceX = locateX1;
		else if (koordSourceX >= sreda[indSreda].x2)
			locateSourceX = locateX2;
		else
			locateSourceX = FindLocate(xNet, nX, koordSourceX);

		if (koordSourceY <= sreda[indSreda].y1)
			locateSourceY = locateY1;
		else if (koordSourceY >= sreda[indSreda].y2)
			locateSourceY = locateY2;
		else
			locateSourceY = FindLocate(yNet, nY, koordSourceY);

		if (koordSourceZ <= sreda[indSreda].z1)
			locateSourceZ = locateZ1;
		else if (koordSourceZ >= sreda[indSreda].z2)
			locateSourceZ = locateZ2;
		else
			locateSourceZ = FindLocate(zNet, nZ, koordSourceZ);

		//CHECK COUNT OF NODES AND INVOKE OPTIMIZATION WITH GREATER NODES
		//RBU
		if (X) OptimizationQuarterX(1, 1, 1, locateSourceX, locateSourceY, locateSourceZ, locateX2, locateY2, locateZ2);
		if (Y) OptimizationQuarterY(1, 1, 1, locateSourceX, locateSourceY, locateSourceZ, locateX2, locateY2, locateZ2);
		if (Z) OptimizationQuarterZ(1, 1, 1, locateSourceX, locateSourceY, locateSourceZ, locateX2, locateY2, locateZ2);
		//LBU
		if (X) OptimizationQuarterX(-1, 1, 1, locateSourceX, locateSourceY, locateSourceZ, locateX1, locateY2, locateZ2);
		if (Y) OptimizationQuarterY(-1, 1, 1, locateSourceX, locateSourceY, locateSourceZ, locateX1, locateY2, locateZ2);
		if (Z) OptimizationQuarterZ(-1, 1, 1, locateSourceX, locateSourceY, locateSourceZ, locateX1, locateY2, locateZ2);
		//RFU
		if (X) OptimizationQuarterX(1, -1, 1, locateSourceX, locateSourceY, locateSourceZ, locateX2, locateY1, locateZ2);
		if (Y) OptimizationQuarterY(1, -1, 1, locateSourceX, locateSourceY, locateSourceZ, locateX2, locateY1, locateZ2);
		if (Z) OptimizationQuarterZ(1, -1, 1, locateSourceX, locateSourceY, locateSourceZ, locateX2, locateY1, locateZ2);
		//LFU
		if (X) OptimizationQuarterX(-1, -1, 1, locateSourceX, locateSourceY, locateSourceZ, locateX1, locateY1, locateZ2);
		if (Y) OptimizationQuarterY(-1, -1, 1, locateSourceX, locateSourceY, locateSourceZ, locateX1, locateY1, locateZ2);
		if (Z) OptimizationQuarterZ(-1, -1, 1, locateSourceX, locateSourceY, locateSourceZ, locateX1, locateY1, locateZ2);

		//RBD
		if (X) OptimizationQuarterX(1, 1, -1, locateSourceX, locateSourceY, locateSourceZ, locateX2, locateY2, locateZ1);
		if (Y) OptimizationQuarterY(1, 1, -1, locateSourceX, locateSourceY, locateSourceZ, locateX2, locateY2, locateZ1);
		if (Z) OptimizationQuarterZ(1, 1, -1, locateSourceX, locateSourceY, locateSourceZ, locateX2, locateY2, locateZ1);
		//LBD
		if (X) OptimizationQuarterX(-1, 1, -1, locateSourceX, locateSourceY, locateSourceZ, locateX1, locateY2, locateZ1);
		if (Y) OptimizationQuarterY(-1, 1, -1, locateSourceX, locateSourceY, locateSourceZ, locateX1, locateY2, locateZ1);
		if (Z) OptimizationQuarterZ(-1, 1, -1, locateSourceX, locateSourceY, locateSourceZ, locateX1, locateY2, locateZ1);
		//RFD
		if (X) OptimizationQuarterX(1, -1, -1, locateSourceX, locateSourceY, locateSourceZ, locateX2, locateY1, locateZ1);
		if (Y) OptimizationQuarterY(1, -1, -1, locateSourceX, locateSourceY, locateSourceZ, locateX2, locateY1, locateZ1);
		if (Z) OptimizationQuarterZ(1, -1, -1, locateSourceX, locateSourceY, locateSourceZ, locateX2, locateY1, locateZ1);
		//LFD
		if (X) OptimizationQuarterX(-1, -1, -1, locateSourceX, locateSourceY, locateSourceZ, locateX1, locateY1, locateZ1);
		if (Y) OptimizationQuarterY(-1, -1, -1, locateSourceX, locateSourceY, locateSourceZ, locateX1, locateY1, locateZ1);
		if (Z) OptimizationQuarterZ(-1, -1, -1, locateSourceX, locateSourceY, locateSourceZ, locateX1, locateY1, locateZ1);
	}
	//duplicatingOfXY();

	profiler << setw(40) << std::left << "Net optimization " <<
		MkaUtils::formattingTime(std::chrono::system_clock::time_point(std::chrono::system_clock::now() - start)) << endl;
}

void Mka3D::PrintLocalMatrix()
{
	return;
	int i, j;
	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < 8; j++)
		{
			logger << setw(15) << localMatrix[i][j];
		}
		logger << endl;
	}
	logger << endl;
}

void Mka3D::PrintPlotMatrix(bool flag_simmeric)
{
	return;
	int i, j;
	double**APlot = new double*[countRegularNodes];
	for (i = 0; i < countRegularNodes; i++)
	{
		APlot[i] = new double[countRegularNodes];
		for (j = 0; j < countRegularNodes; j++)
		{
			APlot[i][j] = 0;
		}
	}
	if (flag_simmeric)
		for (i = 0; i < countRegularNodes; i++)
		{
			APlot[i][i] = di[i];
			for (j = ig[i]; j < ig[i + 1]; j++)
			{
				APlot[i][jg[j]] = ggl[j];
				APlot[jg[j]][i] = ggl[j];
			}
		}
	else
		for (i = 0; i < countRegularNodes; i++)
		{
			APlot[i][i] = di[i];
			for (j = ig[i]; j < ig[i + 1]; j++)
			{
				APlot[i][jg[j]] = ggl[j];
				APlot[jg[j]][i] = ggu[j];
			}
		}

	for (i = 0; i < countRegularNodes; i++)
	{
		for (j = 0; j < countRegularNodes; j++)
		{
			logger << setw(15) << APlot[i][j];
		}
		logger << endl;
	}
	logger << endl;

	for (i = 0; i < countRegularNodes; i++)
	{
		logger << setw(15) << b[i];
	}
	logger << endl;	logger << endl;
}

void Mka3D::inputNet(string sredaInput, string sourceLocate)
{
	int i, numberFields;
	field tempField;
	set<double> xTemp, yTemp, zTemp;
	set<double>::const_iterator it;
	ifstream inpSreda(sredaInput);
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

	ifstream sourceFile(sourceLocate);
	sourceFile >> koordSourceX >> koordSourceY >> koordSourceZ;
	if (koordSourceX > leftX && koordSourceX < rightX)
		xTemp.insert(koordSourceX);
	if (koordSourceY > leftY && koordSourceY < rightY)
		yTemp.insert(koordSourceY);
	if (koordSourceZ > leftZ && koordSourceZ < rightZ)
		zTemp.insert(koordSourceZ);

	if (fabs(areaOfSreda - (rightX - leftX) * (rightY - leftY) * (rightZ - leftZ)) > 1e-15)
	{
		cout << "ќшибка. Ќе правильно заданы подобласти. ќбща€ площадь не равна сумме площадей подобластей." << endl;
		system("pause");
		exit(1);
	}

	nX = xTemp.size();
	delete xNet;
	xNet = new double[nX];
	nY = yTemp.size();
	delete yNet;
	yNet = new double[nY];
	nZ = zTemp.size();
	delete zNet;
	zNet = new double[nZ];
	i = 0;
	for (it = xTemp.begin(); it != xTemp.end(); ++it, i++)
	{
		xNet[i] = *it;
	}
	i = 0;
	for (it = yTemp.begin(); it != yTemp.end(); ++it, i++)
	{
		yNet[i] = *it;
	}
	i = 0;
	for (it = zTemp.begin(); it != zTemp.end(); ++it, i++)
	{
		zNet[i] = *it;
	}
}

void Mka3D::generatePortraitNesoglas()
{
	//int countRegularNodes = xyz_points.size();
	int countLocalIndex = 8;
	set<size_t>* portrait = new set<size_t>[countRegularNodes];
	for (size_t k = 0; k < KE.size(); k++)
	{
		for (size_t i = 0; i < countLocalIndex; i++)
		{
			size_t a = KE[k].uzel[i];
			for (size_t j = 0; j < i; j++)
			{
				size_t b = KE[k].uzel[j];
				// ≈сли оба узла не €вл€ютс€ терминальными
				if (a < countRegularNodes && b < countRegularNodes)
				{
					if (b > a)
						portrait[b].insert(a);
					else
						portrait[a].insert(b);
				}
				else if (a >= countRegularNodes && b < countRegularNodes)
				{
					for (size_t mu = igT[a - countRegularNodes]; mu < igT[a - countRegularNodes + 1]; mu++)
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
				else if (a < countRegularNodes && b >= countRegularNodes)
				{
					for (size_t nu = igT[b - countRegularNodes]; nu < igT[b - countRegularNodes + 1]; nu++)
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
					for (size_t mu = igT[a - countRegularNodes]; mu < igT[a - countRegularNodes + 1]; mu++)
					{
						size_t pos_a = jgT[mu];
						for (size_t nu = igT[b - countRegularNodes]; nu < igT[b - countRegularNodes + 1]; nu++)
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
	for (size_t i = 0; i < countRegularNodes; i++)
		gg_size += portrait[i].size();

	ig = new int[countRegularNodes + 1];
	di = new double[countRegularNodes];
	b = new double[countRegularNodes];
	q = new double[countRegularNodes];
	ig[0] = ig[1] = 0;
	for (size_t i = 0; i < countRegularNodes; i++)
	{
		di[i] = b[i] = q[i] = 0;
		//for (set<size_t>::iterator j = portrait[i].begin(); j != portrait[i].end(); ++j)
		//{
		//	slae.jg[tmp] = *j;
		//	tmp++;
		//}
		ig[i + 1] = ig[i] + portrait[i].size();
	}
	jg = new int[ig[countRegularNodes]];
	ggl = new double[ig[countRegularNodes]];
	ggu = new double[ig[countRegularNodes]];

	for (size_t i = 0; i < ig[countRegularNodes]; i++)
	{
		ggl[i] = 0;
		ggu[i] = 0;
	}

	size_t tmp = 0;
	for (size_t i = 0; i < countRegularNodes; i++)
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
	ofstream igOut(filePrefix + "ig.txt");
	ofstream jgOut(filePrefix + "jg.txt");
	for (i = 0; i <= countRegularNodes; i++)
	{
		igOut << ig[i] << "   " << i << "   ";
		if (i != countRegularNodes)
			igOut << ig[i + 1] - ig[i] << endl;
	}
	cout << endl;
	for (j = 0; j < countRegularNodes; j++)
	{
		jgOut << "string " << j << ": ";
		for (i = ig[j]; i < ig[j + 1]; i++)
		{
			jgOut << jg[i] << " ";
		}
		jgOut << endl;
	}
}

double Mka3D::analiticSolution(Point goal)
{
	return 1 + goal.x + goal.y + goal.z;
	//return goal.x*goal.x + goal.y*goal.y + goal.z*goal.z;
	//return exp(goal.x + goal.y + goal.z);
}

double Mka3D::analiticSolution(double x, double y, double z)
{
	return analiticSolution(Point(x, y, z));
}

double Mka3D::Lambda(int ielem)
{
	return 1;
}

double Mka3D::Gamma(int ielem)
{
	return 1;
}

double Mka3D::Func(int ielem, int node)
{
	return 1 + xyz_points[KE[ielem].uzel[node]].x + xyz_points[KE[ielem].uzel[node]].y + xyz_points[KE[ielem].uzel[node]].z;
	//return xyz_points[KE[ielem].uzel[node]].x*xyz_points[KE[ielem].uzel[node]].x + xyz_points[KE[ielem].uzel[node]].y*xyz_points[KE[ielem].uzel[node]].y +
	//	xyz_points[KE[ielem].uzel[node]].z*xyz_points[KE[ielem].uzel[node]].z - 6;
	//return analiticSolution(xyz_points[KE[ielem].uzel[node]].x, xyz_points[KE[ielem].uzel[node]].y, xyz_points[KE[ielem].uzel[node]].z) - 6;
	//return -2*analiticSolution(xyz_points[KE[ielem].uzel[node]].x, xyz_points[KE[ielem].uzel[node]].y, xyz_points[KE[ielem].uzel[node]].z);
}

void Mka3D::CreateLocalMatrixs(int ielem)
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

void Mka3D::AddToMatrix(int posI, int posJ, double el)
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

void Mka3D::AddToA(int i, int j, double value) {
	if (i < countRegularNodes)
	{
		if (j < countRegularNodes)
		{
			AddToMatrix(i, j, value);
		}
		else
		{
			for (int k = igT[j - countRegularNodes]; k < igT[j - countRegularNodes + 1]; k++)
			{
				AddToMatrix(i, jgT[k], ggT[k] * value);
			}
		}
	}
	else
	{
		if (j < countRegularNodes)
		{
			for (int k = igT[i - countRegularNodes]; k < igT[i - countRegularNodes + 1]; k++)
			{
				AddToMatrix(jgT[k], j, ggT[k] * value);
			}
		}
		else
		{
			for (int k = igT[i - countRegularNodes]; k < igT[i - countRegularNodes + 1]; k++)
			{
				for (int l = igT[j - countRegularNodes]; l < igT[j - countRegularNodes + 1]; l++)
				{
					AddToMatrix(jgT[k], jgT[l], ggT[k] * ggT[l] * value);
				}
			}
		}
	}
}

void Mka3D::AddToB(int i, double value) {
	if (i < countRegularNodes)
	{
		b[i] += value;
	}
	else
		for (int k = igT[i - countRegularNodes]; k < igT[i - countRegularNodes + 1]; k++)
		{
			b[jgT[k]] += value * ggT[k];
		}
}

void Mka3D::Addition(int ielem)
{
	int i, j;
	for (i = 0; i < 8; i++)
	{
		AddToB(KE[ielem].uzel[i], localB[i]);
		for (j = 0; j < 8; j++)
		{
			AddToA(KE[ielem].uzel[i], KE[ielem].uzel[j], localMatrix[i][j]);
		}
	}
}

//intXorYorZ can be 0,1 or 2
void Mka3D::doEdge1(ofstream& outEdge1File, const int intXorYorZ, const int kolvoRegularNode,
	double* varNet1, const int nVarNet1, double* varNet2, const int nVarNet2, const int unknownIndex)
{
	vector<Axis>*info;
	for (int iVar1 = 0; iVar1 < nVarNet1; iVar1++)
	{
		for (int iVar2 = 0; iVar2 < nVarNet2; iVar2++)
		{
			Point *goal;
			bitset<BIT_SIZE> *condition;
			if (intXorYorZ == 0) {
				goal = &Point(xNet[unknownIndex], varNet1[iVar1], varNet2[iVar2]);
				condition = &newNodes[unknownIndex][iVar1][iVar2];
			}
			else if (intXorYorZ == 1) {
				goal = &Point(varNet1[iVar1], yNet[unknownIndex], varNet2[iVar2]);
				condition = &newNodes[iVar1][unknownIndex][iVar2];
			}
			else {
				goal = &Point(varNet1[iVar1], varNet2[iVar2], zNet[unknownIndex]);
				condition = &newNodes[iVar1][iVar2][unknownIndex];
			}
			if (condition->none()) {
				continue;
			}
			int k = indexXYZ(*goal);
			//if (DEBUG) {
				info = nodeInfo(*condition);
			//}
			if (k < kolvoRegularNode) {
				di[k] = 1;
				b[k] = analiticSolution(*goal);
			}
			else {
				continue;
			}
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
			//if (DEBUG) {
				delete info;
			//}
			outEdge1File << k << '\t' << b[k] << endl;
		}
	}
}

void Mka3D::Edge1_not_sim(bool up, bool down, bool left, bool right, bool fore, bool behind)
{
	ofstream ku1(filePrefix + "ku1.txt");
	//int countRegularNodes = xyz_points.size();
	if (down)
		doEdge1(ku1, 2, countRegularNodes, xNet, nX, yNet, nY, 0);
	if (up)
		doEdge1(ku1, 2, countRegularNodes, xNet, nX, yNet, nY, nZ - 1);
	if (left)
		doEdge1(ku1, 0, countRegularNodes, yNet, nY, zNet, nZ, 0);
	if (right)
		doEdge1(ku1, 0, countRegularNodes, yNet, nY, zNet, nZ, nX - 1);
	if (fore)
		doEdge1(ku1, 1, countRegularNodes, xNet, nX, zNet, nZ, 0);
	if (behind)
		doEdge1(ku1, 1, countRegularNodes, xNet, nX, zNet, nZ, nY - 1);
}

int Mka3D::findKE(int ind_nodes[4])
{
	for (size_t indKE = 0; indKE < KE.size(); indKE++)
	{
		int result[4] = { 0,0,0,0 };
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

	cout << "ќшибка findKE. Ќе найден  Ё." << endl;
	system("pause");
	exit(1);
}

//intXorYorZ can be 0,1 or 2
void Mka3D::doEdge2(ofstream& outEdge2File, const int intXorYorZ, const int normalDirect, const int kolvoRegularNode,
	double* varNet1, const int nVarNet1, double* varNet2, const int nVarNet2, const int unknownIndex)
{
	double dUdn[4], h = 0.02, dh1, dh2;
	int indNodes[4];
	for (int iVar1 = 0; iVar1 < nVarNet1; iVar1++)
	{
		for (int iVar2 = 0; iVar2 < nVarNet2; iVar2++)
		{
			int nextVar1, nextVar2;
			Point *dPoint, *gran;
			gran = new Point[4];
			bitset<BIT_SIZE> *condition;
			if (intXorYorZ == 0)
			{
				condition = &newNodes[unknownIndex][iVar1][iVar2];
				if (condition->none() || !canOptimize(-1 * normalDirect, 1, 1, *condition)) {
					continue;
				}
				nextVar1 = findNextY(-1 * normalDirect, 1, 1, iVar1, nVarNet1, unknownIndex, iVar2);
				nextVar2 = findNextZ(-1 * normalDirect, 1, 1, iVar2, nVarNet2, unknownIndex, iVar1);
				if (DEBUG) {
					if (nextVar1 == -1 || nextVar2 == -1
						|| findNextY(-1 * normalDirect, 1, -1, iVar1, nVarNet1, unknownIndex, nextVar2) != nextVar1
						|| findNextZ(-1 * normalDirect, -1, 1, iVar2, nVarNet2, unknownIndex, nextVar1) != nextVar2) {
						logError("Error, no square");
					}
				}
				dPoint = new Point(h, 0, 0);
				gran[0] = Point(xNet[unknownIndex], varNet1[iVar1], varNet2[iVar2]);
				gran[1] = Point(xNet[unknownIndex], varNet1[iVar1], varNet2[nextVar2]);
				gran[2] = Point(xNet[unknownIndex], varNet1[nextVar1], varNet2[iVar2]);
				gran[3] = Point(xNet[unknownIndex], varNet1[nextVar1], varNet2[nextVar2]);
			}
			else if (intXorYorZ == 1) {
				condition = &newNodes[iVar1][unknownIndex][iVar2];
				if (condition->none() || !canOptimize(1, -1 * normalDirect, 1, *condition)) {
					continue;
				}
				nextVar1 = findNextX(1, -1 * normalDirect, 1, iVar1, nVarNet1, unknownIndex, iVar2);
				nextVar2 = findNextZ(1, -1 * normalDirect, 1, iVar2, nVarNet2, iVar1, unknownIndex);
				if (DEBUG) {
					if (nextVar1 == -1 || nextVar2 == -1
						|| findNextX(1, -1 * normalDirect, -1, iVar1, nVarNet1, unknownIndex, nextVar2) != nextVar1
						|| findNextZ(-1, -1 * normalDirect, 1, iVar2, nVarNet2, nextVar1, unknownIndex) != nextVar2) {
						logError("Error, no square");
					}
				}
				dPoint = new Point(0, h, 0);
				gran[0] = Point(varNet1[iVar1], yNet[unknownIndex], varNet2[iVar2]);
				gran[1] = Point(varNet1[iVar1], yNet[unknownIndex], varNet2[nextVar2]);
				gran[2] = Point(varNet1[nextVar1], yNet[unknownIndex], varNet2[iVar2]);
				gran[3] = Point(varNet1[nextVar1], yNet[unknownIndex], varNet2[nextVar2]);
			}
			else
			{
				condition = &newNodes[iVar1][iVar2][unknownIndex];
				if (condition->none() || !canOptimize(1, 1, -1 * normalDirect, *condition)) {
					continue;
				}
				nextVar1 = findNextX(1, 1, -1 * normalDirect, iVar1, nVarNet1, iVar2, unknownIndex);
				nextVar2 = findNextY(1, 1, -1 * normalDirect, iVar2, nVarNet2, iVar1, unknownIndex);
				if (DEBUG) {
					if (nextVar1 == -1 || nextVar2 == -1
						|| findNextX(1, -1, -1 * normalDirect, iVar1, nVarNet1, nextVar2, unknownIndex) != nextVar1
						|| findNextY(-1, 1, -1 * normalDirect, iVar2, nVarNet2, nextVar1, unknownIndex) != nextVar2) {
						logError("Error, no square");
					}
				}
				dPoint = new Point(0, 0, h);
				gran[0] = Point(varNet1[iVar1], varNet2[iVar2], zNet[unknownIndex]);
				gran[1] = Point(varNet1[iVar1], varNet2[nextVar2], zNet[unknownIndex]);
				gran[2] = Point(varNet1[nextVar1], varNet2[iVar2], zNet[unknownIndex]);
				gran[3] = Point(varNet1[nextVar1], varNet2[nextVar2], zNet[unknownIndex]);
			}

			dUdn[0] = normalDirect * (analiticSolution(gran[0] + *dPoint) - analiticSolution(gran[0] - *dPoint)) / (2.0 * h);
			dUdn[1] = normalDirect * (analiticSolution(gran[1] + *dPoint) - analiticSolution(gran[1] - *dPoint)) / (2.0 * h);
			dUdn[2] = normalDirect * (analiticSolution(gran[2] + *dPoint) - analiticSolution(gran[2] - *dPoint)) / (2.0 * h);
			dUdn[3] = normalDirect * (analiticSolution(gran[3] + *dPoint) - analiticSolution(gran[3] - *dPoint)) / (2.0 * h);
			indNodes[0] = indexXYZ(gran[0]);
			indNodes[1] = indexXYZ(gran[1]);
			indNodes[2] = indexXYZ(gran[2]);
			indNodes[3] = indexXYZ(gran[3]);

			int indKE = findKE(indNodes);
			dh1 = varNet1[nextVar1] - varNet1[iVar1];
			dh2 = varNet2[nextVar2] - varNet2[iVar2];
			if (dh1*dh2 == 0) throw new exception("incorrect dh");

			for (size_t linux = 0; linux < 4; linux++)
			{
				double value = 0;
				for (size_t linux2 = 0; linux2 < 4; linux2++)
				{
					value += dUdn[linux2] * M2[linux][linux2] * dh1 * dh2 / 36.0 * Lambda(indKE);
				}
				AddToB(indNodes[linux], value);
				outEdge2File << setw(10) << indNodes[linux] << setw(15) << value << endl;
			}
			delete dPoint;
		}
	}
}

void Mka3D::Edge2_not_sim(bool up, bool down, bool left, bool right, bool fore, bool behind)
{
	ofstream ku2(filePrefix + "ku2.txt");
	//int countRegularNodes = xyz_points.size();
	if (up)
		doEdge2(ku2, 2, 1, countRegularNodes, xNet, nX - 1, yNet, nY - 1, nZ - 1);
	if (down)
		doEdge2(ku2, 2, -1, countRegularNodes, xNet, nX - 1, yNet, nY - 1, 0);
	if (left)
		doEdge2(ku2, 0, -1, countRegularNodes, yNet, nY - 1, zNet, nZ - 1, 0);
	if (right)
		doEdge2(ku2, 0, 1, countRegularNodes, yNet, nY - 1, zNet, nZ - 1, nX - 1);
	if (fore)
		doEdge2(ku2, 1, -1, countRegularNodes, xNet, nX - 1, zNet, nZ - 1, 0);
	if (behind)
		doEdge2(ku2, 1, 1, countRegularNodes, xNet, nX - 1, zNet, nZ - 1, nY - 1);
}

//intXorYorZ can be 0,1 or 2
void Mka3D::doEdge3(ofstream& outEdge3File, int intXorYorZ, int normalDirect, int kolvoRegularNode,
	double* varNet1, int nVarNet1, double* varNet2, int nVarNet2, int unknownIndex)
{
	double dUdn[4], h = 0.01, dh1, dh2, ubetta[4], betta;
	int indNodes[4];
	for (int iVar1 = 0; iVar1 < nVarNet1; iVar1++)
	{
		for (int iVar2 = 0; iVar2 < nVarNet2; iVar2++)
		{
			int nextVar1, nextVar2;
			Point *dPoint, *gran;
			gran = new Point[4];
			bitset<BIT_SIZE> *condition;
			if (intXorYorZ == 0)
			{
				condition = &newNodes[unknownIndex][iVar1][iVar2];
				if (condition->none() || !canOptimize(-1 * normalDirect, 1, 1, *condition)) {
					continue;
				}
				nextVar1 = findNextY(-1 * normalDirect, 1, 1, iVar1, nVarNet1, unknownIndex, iVar2);
				nextVar2 = findNextZ(-1 * normalDirect, 1, 1, iVar2, nVarNet2, unknownIndex, iVar1);
				if (DEBUG) {
					if (nextVar1 == -1 || nextVar2 == -1
						|| findNextY(-1 * normalDirect, 1, -1, iVar1, nVarNet1, unknownIndex, nextVar2) != nextVar1
						|| findNextZ(-1 * normalDirect, -1, 1, iVar2, nVarNet2, unknownIndex, nextVar1) != nextVar2) {
						logError("Error, no square");
					}
				}
				dPoint = new Point(h, 0, 0);
				gran[0] = Point(xNet[unknownIndex], varNet1[iVar1], varNet2[iVar2]);
				gran[1] = Point(xNet[unknownIndex], varNet1[iVar1], varNet2[nextVar2]);
				gran[2] = Point(xNet[unknownIndex], varNet1[nextVar1], varNet2[iVar2]);
				gran[3] = Point(xNet[unknownIndex], varNet1[nextVar1], varNet2[nextVar2]);
			}
			else if (intXorYorZ == 1) {
				condition = &newNodes[iVar1][unknownIndex][iVar2];
				if (condition->none() || !canOptimize(1, -1 * normalDirect, 1, *condition)) {
					continue;
				}
				nextVar1 = findNextX(1, -1 * normalDirect, 1, iVar1, nVarNet1, unknownIndex, iVar2);
				nextVar2 = findNextZ(1, -1 * normalDirect, 1, iVar2, nVarNet2, iVar1, unknownIndex);
				if (DEBUG) {
					if (nextVar1 == -1 || nextVar2 == -1
						|| findNextX(1, -1 * normalDirect, -1, iVar1, nVarNet1, unknownIndex, nextVar2) != nextVar1
						|| findNextZ(-1, -1 * normalDirect, 1, iVar2, nVarNet2, nextVar1, unknownIndex) != nextVar2) {
						logError("Error, no square");
					}
				}
				dPoint = new Point(0, h, 0);
				gran[0] = Point(varNet1[iVar1], yNet[unknownIndex], varNet2[iVar2]);
				gran[1] = Point(varNet1[iVar1], yNet[unknownIndex], varNet2[nextVar2]);
				gran[2] = Point(varNet1[nextVar1], yNet[unknownIndex], varNet2[iVar2]);
				gran[3] = Point(varNet1[nextVar1], yNet[unknownIndex], varNet2[nextVar2]);
			}
			else
			{
				condition = &newNodes[iVar1][iVar2][unknownIndex];
				if (condition->none() || !canOptimize(1, 1, -1 * normalDirect, *condition)) {
					continue;
				}
				nextVar1 = findNextX(1, 1, -1 * normalDirect, iVar1, nVarNet1, iVar2, unknownIndex);
				nextVar2 = findNextY(1, 1, -1 * normalDirect, iVar2, nVarNet2, iVar1, unknownIndex);
				if (DEBUG) {
					if (nextVar1 == -1 || nextVar2 == -1
						|| findNextX(1, -1, -1 * normalDirect, iVar1, nVarNet1, nextVar2, unknownIndex) != nextVar1
						|| findNextY(-1, 1, -1 * normalDirect, iVar2, nVarNet2, nextVar1, unknownIndex) != nextVar2) {
						logError("Error, no square");
					}
				}
				dPoint = new Point(0, 0, h);
				gran[0] = Point(varNet1[iVar1], varNet2[iVar2], zNet[unknownIndex]);
				gran[1] = Point(varNet1[iVar1], varNet2[nextVar2], zNet[unknownIndex]);
				gran[2] = Point(varNet1[nextVar1], varNet2[iVar2], zNet[unknownIndex]);
				gran[3] = Point(varNet1[nextVar1], varNet2[nextVar2], zNet[unknownIndex]);
			}

			dUdn[0] = normalDirect * (analiticSolution(gran[0] + *dPoint) - analiticSolution(gran[0] - *dPoint)) / (2.0 * h);
			dUdn[1] = normalDirect * (analiticSolution(gran[1] + *dPoint) - analiticSolution(gran[1] - *dPoint)) / (2.0 * h);
			dUdn[2] = normalDirect * (analiticSolution(gran[2] + *dPoint) - analiticSolution(gran[2] - *dPoint)) / (2.0 * h);
			dUdn[3] = normalDirect * (analiticSolution(gran[3] + *dPoint) - analiticSolution(gran[3] - *dPoint)) / (2.0 * h);
			indNodes[0] = indexXYZ(gran[0]);
			indNodes[1] = indexXYZ(gran[1]);
			indNodes[2] = indexXYZ(gran[2]);
			indNodes[3] = indexXYZ(gran[3]);

			int indKE = findKE(indNodes);
			dh1 = varNet1[nextVar1] - varNet1[iVar1];
			dh2 = varNet2[nextVar2] - varNet2[iVar2];
			if (dh1*dh2 == 0) throw new exception("incorrect dh");
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
					AddToA(indNodes[linux], indNodes[linux2], betta * M2[linux][linux2] * dh1 * dh2 / 36.0);
				}
				AddToB(indNodes[linux], betta * value);
				outEdge3File << setw(10) << indNodes[linux] << setw(15) << value << endl;
			}
		}
	}
}

void Mka3D::Edge3_not_sim(bool up, bool down, bool left, bool right, bool fore, bool behind)
{
	ofstream ku3(filePrefix + "ku3.txt");
	//int countRegularNodes = xyz_points.size();
	if (up)
		doEdge3(ku3, 2, 1, countRegularNodes, xNet, nX - 1, yNet, nY - 1, nZ - 1);
	if (down)
		doEdge3(ku3, 2, -1, countRegularNodes, xNet, nX - 1, yNet, nY - 1, 0);
	if (left)
		doEdge3(ku3, 0, -1, countRegularNodes, yNet, nY - 1, zNet, nZ - 1, 0);
	if (right)
		doEdge3(ku3, 0, 1, countRegularNodes, yNet, nY - 1, zNet, nZ - 1, nX - 1);
	if (fore)
		doEdge3(ku3, 1, -1, countRegularNodes, xNet, nX - 1, zNet, nZ - 1, 0);
	if (behind)
		doEdge3(ku3, 1, 1, countRegularNodes, xNet, nX - 1, zNet, nZ - 1, nY - 1);
}

void Mka3D::generateGlobalMatrix()
{
	auto start = std::chrono::system_clock::now();

	int ielem, i;
	//int countRegularNodes = xyz_points.size();
	for (ielem = 0; ielem < KE.size(); ielem++)
	{
		CreateLocalMatrixs(ielem);
		Addition(ielem);
		PrintLocalMatrix();
	}
	//PrintPlotMatrix(false);

	Edge2_not_sim(1, 1, 1, 1, 1, 1);
	//Edge3_not_sim(1, 1, 1, 1, 1, 1);
	for (i = 0; i < ig[countRegularNodes]; i++)
	{
		ggu[i] = ggl[i];
	}
	//Edge1_not_sim(1, 1, 0, 0, 0, 0);
	//Edge1_not_sim(1, 1, 1, 1, 1, 1);

	profiler << setw(40) << std::left << "Generate global matrix " <<
		MkaUtils::formattingTime(std::chrono::system_clock::time_point(std::chrono::system_clock::now() - start)) << endl;
}

void Mka3D::mult(double* res, double* v)
{
	//int countRegularNodes = xyz_points.size();
	for (int i = 0; i < countRegularNodes; i++)
		res[i] = 0;
	for (int i = 0; i < countRegularNodes; i++)
	{
		for (int j = ig[i]; j < ig[i + 1]; j++)
		{
			res[i] += ggl[j] * v[jg[j]];
			res[jg[j]] += ggu[j] * v[i];
		}
		res[i] += di[i] * v[i];
	}
}

double Mka3D::ScalarMult(double* v1, double* v2)
{
	int i;
	double result;
	//int countRegularNodes = xyz_points.size();
	result = 0;
	for (i = 0; i < countRegularNodes; i++)
	{
		result += v1[i] * v2[i];
	}
	return result;
}

//	”множение матрицы на вектор
void Mka3D::MultMatrixOnVector(double* in, double* out)
{
	int i, j;
	double* out1;
	//int countRegularNodes = xyz_points.size();
	out1 = new double[countRegularNodes];
	for (i = 0; i < countRegularNodes; i++)
	{
		out1[i] = di[i] * in[i];
		for (j = ig[i]; j < ig[i + 1]; j++)
		{
			out1[i] += ggl[j] * in[jg[j]];
			out1[jg[j]] += ggu[j] * in[i];
		}
	}
	for (i = 0; i < countRegularNodes; i++)
		out[i] = out1[i];
	delete[] out1;
}

void Mka3D::calcPogreshnost(ofstream& output)
{
	int i;
	double sumCh = 0, sumZn = 0;
	output << setw(10) << "x" << setw(10) << "y" << setw(10) << "z" << setw(18) << "analitic" << setw(18) << "solution" << setw(18) << "pogreshn" << endl;
	for (i = 0; i < countRegularNodes; i++)
	{
		sumCh += (q[i] - analiticSolution(xyz_points[i])) * (q[i] - analiticSolution(xyz_points[i]));
		output << setw(10) << xyz_points[i].x << setw(10) << xyz_points[i].y << setw(10) << xyz_points[i].z << setw(18) << analiticSolution(xyz_points[i]) << setw(18) << q[i] << setw(18) << q[i] - analiticSolution(xyz_points[i]) << endl;
		sumZn += analiticSolution(xyz_points[i]) * analiticSolution(xyz_points[i]);
	}
	output << "KE number = " << KE.size() << endl;
	output << "Nodes number = " << countRegularNodes << endl;
	output << endl << "ќтносительна€ погрешность  = " << sqrt(sumCh / sumZn);
}

void Mka3D::runLOS()
{
	auto start = std::chrono::system_clock::now();

	int maxiter = 10000, i;
	double alfa, alfachisl, alfaznam, beta, betachisl, betaznam, checkE, epsMSG = 1e-16;

	//int countRegularNodes = xyz_points.size();
	double* r = new double[countRegularNodes];
	double* s = new double[countRegularNodes];
	double* z = new double[countRegularNodes];
	double* p = new double[countRegularNodes];
	double* rout = new double[countRegularNodes];

	for (i = 0; i < countRegularNodes; i++)
	{
		s[i] = rout[i] = r[i] = q[i] = z[i] = p[i] = 0;
	}
	MultMatrixOnVector(q, r);
	for (i = 0; i < countRegularNodes; i++)
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
		for (i = 0; i < countRegularNodes; i++)
		{
			q[i] = q[i] + alfa * z[i];
			r[i] = r[i] - alfa * p[i];
		}
		MultMatrixOnVector(r, rout);
		betachisl = ScalarMult(p, rout);
		betaznam = ScalarMult(p, p);
		beta = -betachisl / betaznam;
		for (i = 0; i < countRegularNodes; i++)
		{
			z[i] = r[i] + beta * z[i];
			p[i] = rout[i] + beta * p[i];
		}
		checkE = sqrt(ScalarMult(r, r) / ScalarMult(b, b));
	}
	profiler << setw(40) << std::left << "Solve slau " <<
		MkaUtils::formattingTime(std::chrono::system_clock::time_point(std::chrono::system_clock::now() - start)) << endl;
}

Mka3D::locateOfPoint Mka3D::FindLocate(Point sample)
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
		throw new exception("ќшибка FindLocate: не найдена точка");
	}
	return a;
}

void Mka3D::sigmTChain(int nTermNode, int startOfChain, double mnojT, set<int> &visitedNodes)
{
	int i, j;
	for each (neighbor var in tmpSigm[startOfChain].neighbors)
	{
		//if (!WithoutPerehlest)
		//for (i = 0; i < sigmT[nTermNode].neighbors.size(); i++)
		//{
		//	if (sigmT[nTermNode].neighbors[i] == sigmT[startOfChain].neighbors[j])
		//	{
		//		cout << "ќбнаружена цепочка терминальных узлов в " << nTermNode << endl;
		//		//system("pause");
		//		//exit(1);
		//	}
		//}
		if (visitedNodes.find(var.index) != visitedNodes.end()) {
			logger << "found circle " << var.index << " in " << startOfChain << endl;
			continue;
		}
		visitedNodes.insert(var.index);
		if (var.index >= countRegularNodes)
		{
			sigmTChain(nTermNode, var.index - countRegularNodes, mnojT * var.weight, visitedNodes);
		}
		else {
			sigmNewT[nTermNode].neighbors.insert(neighbor(var.index, mnojT*var.weight));
		}
	}
}

set<byte> Mka3D::getValues(pair<map<Point, byte>::iterator, map<Point, byte>::iterator> &range) {
	set<byte> result;
	for (map<Point, byte>::iterator it = range.first; it != range.second; it++)
	{
		byte e = it->second;
		result.insert(e);
	}
	return result;
}

void Mka3D::genT3D() {
	ofstream outT(filePrefix + "Tmatrix.txt");
	igT = new int[nColT + 1];

	//формируем вспомогательную сигм-структуру
	//sigmNewT = new sigmStruct3D[nColT];
	sigmNewT.resize(nColT);
	tmpSigm = new sigmStruct3D[nColT];

	for (int i = countRegularNodes; i < xyz_points.size(); i++)
	{
		int inc = i - countRegularNodes;
		tmpSigm[inc].terminalNode = i;
		Point target = xyz_points[i];
		locateOfPoint term = FindLocate(xyz_points[i]);
		if (newNodes[term.i][term.j][term.k].test(IS_REGULAR))
			throw new exception("wrong terminal node in genT3D");

		pair<map<Point, byte>::iterator, map<Point, byte>::iterator> range
			= termNodeOnEdge.equal_range(target);
		set<byte> edges = getValues(range);
		
		if (edges.size() > 1) {
			logError("Hey dude, i cant calc T with few term nodes in one line");
		}
		if (edges.find(X_EDGE_DELETED) != edges.end())
		{
			int neibLeft;
			for (neibLeft = term.i - 1; neibLeft >= 0; neibLeft--)
			{
				if (newNodes[neibLeft][term.j][term.k].test(IS_REGULAR) ||
					(!newNodes[neibLeft][term.j][term.k].none() && !hasLeft(newNodes[neibLeft][term.j][term.k]) ||
						getValues(termNodeOnEdge.equal_range(Point(xNet[neibLeft], yNet[term.j], zNet[term.k]))).count(X_EDGE_DELETED) == 0)) 
					break;
			}
			int neibRight;
			for (neibRight = term.i + 1; neibRight < nX; neibRight++)
			{
				if (newNodes[neibRight][term.j][term.k].test(IS_REGULAR) ||
					(!newNodes[neibRight][term.j][term.k].none() && !hasRight(newNodes[neibRight][term.j][term.k]) ||
						getValues(termNodeOnEdge.equal_range(Point(xNet[neibRight], yNet[term.j], zNet[term.k]))).count(X_EDGE_DELETED) == 0))
					break;
			}
			double length = xNet[neibRight] - xNet[neibLeft];
			int iL = indexXYZ(xNet[neibLeft], yNet[term.j], zNet[term.k]);
			int iR = indexXYZ(xNet[neibRight], yNet[term.j], zNet[term.k]);
			tmpSigm[inc].neighbors.insert(neighbor(iL, (xNet[neibRight] - xNet[term.i]) / length));
			tmpSigm[inc].neighbors.insert(neighbor(iR, (xNet[term.i] - xNet[neibLeft]) / length));
		}
		if (edges.find(Y_EDGE_DELETED) != edges.end())
		{
			int neibBack;
			for (neibBack = term.j + 1; neibBack < nY; neibBack++)
			{
				if (newNodes[term.i][neibBack][term.k].test(IS_REGULAR) ||
					(!newNodes[term.i][neibBack][term.k].none() && !hasBack(newNodes[term.i][neibBack][term.k]) ||
						getValues(termNodeOnEdge.equal_range(Point(xNet[term.i], yNet[neibBack], zNet[term.k]))).count(Y_EDGE_DELETED) == 0))
					break;
			}
			int neibFore;
			for (neibFore = term.j - 1; neibFore >= 0; neibFore--)
			{
				if (newNodes[term.i][neibFore][term.k].test(IS_REGULAR) ||
					(!newNodes[term.i][neibFore][term.k].none() && !hasFore(newNodes[term.i][neibFore][term.k]) ||
						getValues(termNodeOnEdge.equal_range(Point(xNet[term.i], yNet[neibFore], zNet[term.k]))).count(Y_EDGE_DELETED) == 0)) 
					break;
			}
			double length = yNet[neibBack] - yNet[neibFore];
			int iF = indexXYZ(xNet[term.i], yNet[neibFore], zNet[term.k]);
			int iB = indexXYZ(xNet[term.i], yNet[neibBack], zNet[term.k]);
			tmpSigm[inc].neighbors.insert(neighbor(iF, (yNet[neibBack] - yNet[term.j]) / length));
			tmpSigm[inc].neighbors.insert(neighbor(iB, (yNet[term.j] - yNet[neibFore]) / length));
		}
		if (edges.find(Z_EDGE_DELETED) != edges.end())
		{
			int neibUp;
			for (neibUp = term.k + 1; neibUp < nZ; neibUp++)
			{
				if (newNodes[term.i][term.j][neibUp].test(IS_REGULAR) ||
					(!newNodes[term.i][term.j][neibUp].none() && !hasUp(newNodes[term.i][term.j][neibUp]) ||
						getValues(termNodeOnEdge.equal_range(Point(xNet[term.i], yNet[term.j], zNet[neibUp]))).count(Z_EDGE_DELETED) == 0)) 
					break;
			}
			int neibDown;
			for (neibDown = term.k - 1; neibDown >= 0; neibDown--)
			{
				if (newNodes[term.i][term.j][neibDown].test(IS_REGULAR) ||
					(!newNodes[term.i][term.j][neibDown].none() && !hasDown(newNodes[term.i][term.j][neibDown]) ||
						getValues(termNodeOnEdge.equal_range(Point(xNet[term.i], yNet[term.j], zNet[neibDown]))).count(Z_EDGE_DELETED) == 0)) 
					break;
			}
			double length = zNet[neibUp] - zNet[neibDown];
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
		for each (neighbor var in tmpSigm[i].neighbors)
		{
			visitedNodes.insert(var.index);
			if (var.index >= countRegularNodes) {
				sigmTChain(i, var.index - countRegularNodes, var.weight, visitedNodes);
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

bool Mka3D::belongToInterval(double val, double l1, double l2) {
	return MkaUtils::compare(val, l1) >= 0 && MkaUtils::compare(val, l2) <= 0;
}

bool Mka3D::containsPointOnVergeOfKe(int termNode, int iKe, vector<pair<Point, Point>>* pairs) {
	Point& target = xyz_points[termNode];
	bool on_x_line = MkaUtils::equals(target.x, xyz_points[KE[iKe].uzel[0]].x) || MkaUtils::equals(target.x, xyz_points[KE[iKe].uzel[1]].x);
	bool on_y_line = MkaUtils::equals(target.y, xyz_points[KE[iKe].uzel[0]].y) || MkaUtils::equals(target.y, xyz_points[KE[iKe].uzel[2]].y);
	bool on_z_line = MkaUtils::equals(target.z, xyz_points[KE[iKe].uzel[0]].z) || MkaUtils::equals(target.z, xyz_points[KE[iKe].uzel[4]].z);
	if (!(on_x_line || on_y_line || on_z_line)) return false;
	if (on_x_line && on_y_line && on_z_line) return false;
	vector<Axis>*info;
	for (int i = 0; i < 8; i++)
	{
		if (xyz_points[KE[iKe].uzel[i]] == target) {
			return false;
		}
	}
	bool error = false;
	//do {
		int inc = termNode - countRegularNodes;
		tmpSigm[inc].terminalNode = termNode;
		if (on_x_line) {
			if (belongToInterval(target.y, xyz_points[KE[iKe].uzel[0]].y, xyz_points[KE[iKe].uzel[2]].y) &&
				belongToInterval(target.z, xyz_points[KE[iKe].uzel[0]].z, xyz_points[KE[iKe].uzel[4]].z)) {
				locateOfPoint loc = FindLocate(target);
				info = nodeInfo(newNodes[loc.i][loc.j][loc.k]);
				if (hasYLines(newNodes[loc.i][loc.j][loc.k])) {
					int i1 = indexXYZ(target.x, xyz_points[KE[iKe].uzel[0]].y, target.z);
					int i2 = indexXYZ(target.x, xyz_points[KE[iKe].uzel[2]].y, target.z);
					if (i1 != -1 && i2 != -1) {
						double length = xyz_points[KE[iKe].uzel[2]].y - xyz_points[KE[iKe].uzel[0]].y;
						pairs[inc].push_back(pair<Point, Point>(
							Point(target.x, xyz_points[KE[iKe].uzel[0]].y, target.z),
							Point(target.x, xyz_points[KE[iKe].uzel[2]].y, target.z))
						);
						tmpSigm[inc].neighbors.insert(neighbor(i1, (xyz_points[KE[iKe].uzel[2]].y - target.y) / length));
						tmpSigm[inc].neighbors.insert(neighbor(i2, (target.y - xyz_points[KE[iKe].uzel[0]].y) / length));
						return true;
					}
				}
				if (hasZLines(newNodes[loc.i][loc.j][loc.k])) {
					int i1 = indexXYZ(target.x, target.y, xyz_points[KE[iKe].uzel[0]].z);
					int i2 = indexXYZ(target.x, target.y, xyz_points[KE[iKe].uzel[4]].z);
					if (i1 != -1 && i2 != -1) {
						double length = xyz_points[KE[iKe].uzel[4]].z - xyz_points[KE[iKe].uzel[0]].z;
						pairs[inc].push_back(pair<Point, Point>(
							Point(target.x, target.y, xyz_points[KE[iKe].uzel[0]].z),
							Point(target.x, target.y, xyz_points[KE[iKe].uzel[4]].z))
						);
						tmpSigm[inc].neighbors.insert(neighbor(i1, (xyz_points[KE[iKe].uzel[2]].z - target.z) / length));
						tmpSigm[inc].neighbors.insert(neighbor(i2, (target.z - xyz_points[KE[iKe].uzel[0]].z) / length));
						return true;
					}
				}
				error = true;
				delete info;
				logError("WTF");
			}
		}
		if (on_y_line) {
			if (belongToInterval(target.x, xyz_points[KE[iKe].uzel[0]].x, xyz_points[KE[iKe].uzel[1]].x) &&
				belongToInterval(target.z, xyz_points[KE[iKe].uzel[0]].z, xyz_points[KE[iKe].uzel[4]].z)) {
				locateOfPoint loc = FindLocate(target);
				info = nodeInfo(newNodes[loc.i][loc.j][loc.k]);
				if (hasXLines(newNodes[loc.i][loc.j][loc.k])) {
					int i1 = indexXYZ(xyz_points[KE[iKe].uzel[0]].x, target.y, target.z);
					int i2 = indexXYZ(xyz_points[KE[iKe].uzel[1]].x, target.y, target.z);
					if (i1 != -1 && i2 != -1) {
						double length = xyz_points[KE[iKe].uzel[1]].x - xyz_points[KE[iKe].uzel[0]].x;
						pairs[inc].push_back(pair<Point, Point>(
							Point(xyz_points[KE[iKe].uzel[0]].x, target.y, target.z),
							Point(xyz_points[KE[iKe].uzel[1]].x, target.y, target.z))
						);
						tmpSigm[inc].neighbors.insert(neighbor(i1, (xyz_points[KE[iKe].uzel[1]].x - target.x) / length));
						tmpSigm[inc].neighbors.insert(neighbor(i2, (target.x - xyz_points[KE[iKe].uzel[0]].x) / length));
						return true;
					}
				}
				if (hasZLines(newNodes[loc.i][loc.j][loc.k])) {
					int i1 = indexXYZ(target.x, target.y, xyz_points[KE[iKe].uzel[0]].z);
					int i2 = indexXYZ(target.x, target.y, xyz_points[KE[iKe].uzel[4]].z);
					if (i1 != -1 && i2 != -1) {
						double length = xyz_points[KE[iKe].uzel[4]].z - xyz_points[KE[iKe].uzel[0]].z;
						pairs[inc].push_back(pair<Point, Point>(
							Point(target.x, target.y, xyz_points[KE[iKe].uzel[0]].z),
							Point(target.x, target.y, xyz_points[KE[iKe].uzel[4]].z))
						);
						tmpSigm[inc].neighbors.insert(neighbor(i1, (xyz_points[KE[iKe].uzel[2]].z - target.z) / length));
						tmpSigm[inc].neighbors.insert(neighbor(i2, (target.z - xyz_points[KE[iKe].uzel[0]].z) / length));
						return true;
					}
				}
				error = true;
				delete info;
				logError("WTF");
			}
		}
		if (on_z_line) {
			if (belongToInterval(target.x, xyz_points[KE[iKe].uzel[0]].x, xyz_points[KE[iKe].uzel[1]].x) &&
				belongToInterval(target.y, xyz_points[KE[iKe].uzel[0]].y, xyz_points[KE[iKe].uzel[2]].y)) {
				locateOfPoint loc = FindLocate(target);
				info = nodeInfo(newNodes[loc.i][loc.j][loc.k]);
				if (hasXLines(newNodes[loc.i][loc.j][loc.k])) {
					int i1 = indexXYZ(xyz_points[KE[iKe].uzel[0]].x, target.y, target.z);
					int i2 = indexXYZ(xyz_points[KE[iKe].uzel[1]].x, target.y, target.z);
					if (i1 != -1 && i2 != -1) {
						double length = xyz_points[KE[iKe].uzel[1]].x - xyz_points[KE[iKe].uzel[0]].x;
						pairs[inc].push_back(pair<Point, Point>(
							Point(xyz_points[KE[iKe].uzel[0]].x, target.y, target.z),
							Point(xyz_points[KE[iKe].uzel[1]].x, target.y, target.z))
						);
						tmpSigm[inc].neighbors.insert(neighbor(i1, (xyz_points[KE[iKe].uzel[1]].x - target.x) / length));
						tmpSigm[inc].neighbors.insert(neighbor(i2, (target.x - xyz_points[KE[iKe].uzel[0]].x) / length));
						return true;
					}
				}
				if (hasYLines(newNodes[loc.i][loc.j][loc.k])) {
					int i1 = indexXYZ(target.x, xyz_points[KE[iKe].uzel[0]].y, target.z);
					int i2 = indexXYZ(target.x, xyz_points[KE[iKe].uzel[2]].y, target.z);
					if (i1 != -1 && i2 != -1) {
						double length = xyz_points[KE[iKe].uzel[2]].y - xyz_points[KE[iKe].uzel[0]].y;
						pairs[inc].push_back(pair<Point, Point>(
							Point(target.x, xyz_points[KE[iKe].uzel[0]].y, target.z),
							Point(target.x, xyz_points[KE[iKe].uzel[2]].y, target.z))
						);
						tmpSigm[inc].neighbors.insert(neighbor(i1, (xyz_points[KE[iKe].uzel[2]].y - target.y) / length));
						tmpSigm[inc].neighbors.insert(neighbor(i2, (target.y - xyz_points[KE[iKe].uzel[0]].y) / length));
						return true;
					}
				}
				error = true;
				delete info;
				logError("WTF");
			}
		}
	return false;
}

void Mka3D::constructXyzAndNvtr()
{
	int i, j, t, k;
	ofstream fileXY(filePrefix + "xyz.txt");
	ofstream fileNvtr(filePrefix + "nvtr.txt");
	vector<Point> ncPoint;

	//формируем массивы регул€рных и терминальных вершин
	for (t = 0; t < nZ; t++)
	{
		for (j = 0; j < nY; j++)
		{
			for (i = 0; i < nX; i++)
			{
				vector<Axis>* info;
				//if (DEBUG) {
					info = nodeInfo(newNodes[i][j][t]);
					//if (!newNodes[i][j][t].test(IS_REGULAR) && hasAll(newNodes[i][j][t]))
					//	throw new exception("it is happened!");
				//}
				if (newNodes[i][j][t].test(IS_REGULAR)/* || hasAll(newNodes[i][j][t])*/)	//может быть и терминальным, добавить проверку
				{
					xyz_points.push_back(Point(xNet[i], yNet[j], zNet[t]));
				}
				else if (!newNodes[i][j][t].none())
				{
					ncPoint.push_back(Point(xNet[i], yNet[j], zNet[t]));
				}
				//if (DEBUG) {
					delete info;
				//}
			}
		}
	}

	//добавл€ем в конец массива вершин терминальные вершины
	for (vector<Point>::iterator it = ncPoint.begin(); it < ncPoint.end(); it++)
	{
		xyz_points.push_back(*it);
	}

	//количество столбцов в T
	nColT = ncPoint.size();
	countRegularNodes = xyz_points.size() - nColT;

	//формируем файл xyz.txt
	fileXY << xyz_points.size() << " " << countRegularNodes << endl;
	for (i = 0; i < xyz_points.size(); i++)
	{
		fileXY << xyz_points[i].x << " " << xyz_points[i].y << " " << xyz_points[i].z << endl;
	}

	//формируем структуру  Ё
	nvtr tempNvtr;
	int countKe = 0;
	//tmpSigm = new sigmStruct3D[nColT];
	//vector<pair<Point, Point>>* vectorPairs = new vector<pair<Point, Point>>[nColT];
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

					tempNvtr.numberField = FindAreaNumber(tempNvtr.uzel);
					KE.push_back(tempNvtr);

					//for (int i = countRegularNodes; i < xyz_points.size(); i++)
					//{
					//	containsPointOnVergeOfKe(i, countKe, vectorPairs);
					//}
					countKe++;
				}
			}
		}
	}

	//for (int i = 0; i < nColT; i++)
	//{
	//	set<int> visitedNodes;
	//	visitedNodes.insert(tmpSigm[i].terminalNode);
	//	sigmNewT[i].terminalNode = tmpSigm[i].terminalNode;
	//	logger << "Start chaining " << tmpSigm[i].terminalNode << " node" << endl;
	//	for each (neighbor var in tmpSigm[i].neighbors)
	//	{
	//		visitedNodes.insert(var.index);
	//		if (var.index >= countRegularNodes) {
	//			sigmTChain(i, var.index - countRegularNodes, var.weight, visitedNodes);
	//		}
	//		else {
	//			sigmNewT[i].neighbors.insert(var);
	//		}
	//	}
	//	//check
	//	double sum = 0;
	//	for each (neighbor var in sigmNewT[i].neighbors)
	//	{
	//		sum += var.weight;
	//	}
	//	if (!MkaUtils::equals(sum, 1)) {
	//		logError("Incorrect chain of matrix T, weights sum not equals to 1");
	//	}
	//}

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

void Mka3D::prepareNetForSolve() {

}
