#pragma once
#include "MkaUtils.h"

using namespace std;

class Point_cylindrical
{
public:
	double r, z;
	Point_cylindrical();
	Point_cylindrical(double r,double z);
	~Point_cylindrical();
	friend bool operator<(const Point_cylindrical & p_lhs, const Point_cylindrical & p_rhs);
	friend ostream& operator<<(ostream& os, const Point_cylindrical& point);
	friend bool operator==(const Point_cylindrical& lhs, const Point_cylindrical& rhs);
	friend bool operator!=(const Point_cylindrical& lhs, const Point_cylindrical& rhs);
	Point_cylindrical operator- (Point_cylindrical& left);
	Point_cylindrical operator+ (Point_cylindrical& left);

private:

};

