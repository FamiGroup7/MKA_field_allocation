#pragma once
#include "Point.h"

class Edge
{
public:
	Edge();
	Edge(Point p1, Point p2);
	~Edge();
	friend bool operator<(const Edge & one, const Edge & two);
	friend bool operator==(const Edge& lhs, const Edge& rhs);
	friend bool operator!=(const Edge& lhs, const Edge& rhs);
	friend ostream& operator<<(ostream& os, const Edge& edge);
	Point p1, p2;
};
