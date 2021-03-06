#include "Point.h"

Point::Point()
{
}

Point::Point(double x, double y, double z) : Point(x, y, z, -1)
{
}

Point::Point(double x, double y, double z, int i) : x(x), y(y), z(z), ind(i)
{
}

Point::~Point()
{
}

bool operator < (const Point & left, const Point & right)
{
	int resZ = MkaUtils::compare(left.z, right.z);
	if (resZ != 0)return resZ < 0;
	int resY = MkaUtils::compare(left.y, right.y);
	if (resY != 0)return resY < 0;
	int resX = MkaUtils::compare(left.x, right.x);
	return resX < 0;
	//return MkaUtils::compare(left.z, right.z) <= 0
	//	&& MkaUtils::compare(left.y, right.y) <= 0
	//	&& MkaUtils::compare(left.x, right.x) <= 0;
}

Point Point::operator-(Point & right)
{
	return Point(this->x - right.x, this->y - right.y, this->z - right.z);
}

Point Point::operator+(Point & right)
{
	return Point(right.x + this->x, right.y + this->y, right.z + this->z);
}

double & Point::operator[](size_t i)
{
	if (i == 0) return this->x;
	if (i == 1) return this->y;
	if (i == 2) return this->z;
	throw new exception("Unavailable index of point");
}

ostream & operator<<(ostream & os, const Point & point)
{
	os << "{ x = " << point.x << "; y = " << point.y << "; z = " << point.z << "}";
	return os;
}

bool operator==(const Point& lhs, const Point& rhs)
{
	return MkaUtils::equals(lhs.x, rhs.x)
		&& MkaUtils::equals(lhs.y, rhs.y)
		&& MkaUtils::equals(lhs.z, rhs.z);
}

bool operator!=(const Point& lhs, const Point& rhs)
{
	return !(lhs == rhs);
}

bool Point::inside(double x0, double x1, double y0, double y1, double z0, double z1) const
{
	return x >= x0 && x <= x1 && y >= y0 && y <= y1 && z >= z0 && z <= z1;
}