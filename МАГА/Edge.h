#pragma once
#include "Point.h"

class Edge
{
public:
	Edge();
	Edge(Point p1, Point p2);
	~Edge();
	Point p1, p2;
};
