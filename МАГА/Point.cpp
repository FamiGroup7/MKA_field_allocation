#include "Point.h"

Point::Point()
{
}

Point::Point(double x, double y, double z) :x(x), y(y), z(z)
{
}

Point::~Point()
{
}

bool operator < (const Point & left, const Point & right)
{
	return MkaUtils::compare(left.z, right.z) < 0
		|| MkaUtils::compare(left.y, right.y) < 0
		|| MkaUtils::compare(left.x, right.x) < 0;
}