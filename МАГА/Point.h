#pragma once
#include "MkaUtils.h"

class Point 
{
public:
	Point();
	Point(double x, double y, double z);
	~Point();
	double x, y, z;
	friend bool operator<(const Point & p_lhs, const Point & p_rhs);
private:

};