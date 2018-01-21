#pragma once
#include "MkaUtils.h"
#include "Point.h"

#define HELP_STRUCT_KE_MAX_SIZE 50

struct field
{
	double x1, x2, y1, y2, z1, z2, lambda, gamma;

	field(double x1, double x2, double y1, double y2, double z1, double z2, double lambda, double gamma) :
		x1(x1), x2(x2), y1(y1), y2(y2), z1(z1), z2(z2),
		lambda(lambda), gamma(gamma)
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
		memcpy(uzel, newNodes, sizeof(int) * 8);
	}

	nvtr() {}
};

struct locateOfPoint
{
	int i, j, k;

	locateOfPoint(int i, int j, int k) : i(i), j(j), k(k) {}

	locateOfPoint() : i(0), j(0), k(0) {}
};

class ke_tree {
public:

	ke_tree();
	~ke_tree();
	void init(double x1, double x2, double y1, double y2, double z1, double z2, vector<Point>*points, vector<nvtr>* kes);
	Point centerOfKe(int iKe);
	void addKe(int iKe);
	int findKe(Point point);

private:
	int delim;
	double x1_outer, x2_outer, y1_outer, y2_outer, z1_outer, z2_outer;
	vector<int> elements;
	ke_tree*child1 = NULL, *child2 = NULL;
	field info;
	vector<Point>* points;
	vector<nvtr>* kes;

	void expandOuter(double probablyNewOuter, double&outer, int comparison);
	void init(int parentDelim, ke_tree &parent, bool first);
	bool contains(Point point);
	bool containsInOuter(Point point);
	bool addKe(int iKe, Point point);
	void expand(Point leftDown, Point rightUp);
};