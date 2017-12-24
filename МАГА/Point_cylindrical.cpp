#include "Point_cylindrical.h"

Point_cylindrical::Point_cylindrical()
{
}

Point_cylindrical::Point_cylindrical(double r, double z) :r(r), z(z)
{
}

Point_cylindrical::~Point_cylindrical()
{
}

bool operator < (const Point_cylindrical & left, const Point_cylindrical & right)
{
	int resY = MkaUtils::compare(left.z, right.z);
	if (resY != 0)return resY < 0;
	int resX = MkaUtils::compare(left.r, right.r);
	return resX < 0;
}

Point_cylindrical Point_cylindrical::operator-(Point_cylindrical & right)
{
	return Point_cylindrical(this->r - right.r, this->z - right.z);
}

Point_cylindrical Point_cylindrical::operator+(Point_cylindrical & right)
{
	return Point_cylindrical(right.r + this->r, right.z + this->z);
}

ostream & operator<<(ostream & os, const Point_cylindrical & point)
{
	os << "{ r = " << point.r << "; z = " << point.z << "}";
	return os;
}

bool operator==(const Point_cylindrical& lhs, const Point_cylindrical& rhs)
{
	return MkaUtils::equals(lhs.r, rhs.r)
		&& MkaUtils::equals(lhs.z, rhs.z);
}

bool operator!=(const Point_cylindrical& lhs, const Point_cylindrical& rhs)
{
	return !(lhs == rhs);
}
