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