#include "Edge.h"

Edge::Edge()
{
}

Edge::Edge(Point p1, Point p2) {
	if (p1 < p2) {
		this->p1 = p1;
		this->p2 = p2;
	}
	else {
		this->p1 = p2;
		this->p2 = p1;
	}
}

Edge::~Edge()
{
}

bool operator<(const Edge & one, const Edge & two) {
	if (one.p1 < two.p1) return true;
	else if (one.p1 != two.p1)return false;
	return one.p2 < two.p2;
}

bool operator==(const Edge& lhs, const Edge& rhs)
{
	return lhs.p1 == rhs.p1 && lhs.p2 == rhs.p2
		|| lhs.p1 == rhs.p2 && lhs.p2 == rhs.p1;
}

bool operator!=(const Edge& lhs, const Edge& rhs)
{
	return !(lhs == rhs);
}

ostream & operator<<(ostream & os, const Edge & edge)
{
	os << "{\n\tp1 = " << edge.p1 << "\n\tp2 = " << edge.p2 << "\n}\n";
	return os;
}