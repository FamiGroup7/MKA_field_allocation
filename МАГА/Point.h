#pragma once
#include "MkaUtils.h"
#include <iostream>

using namespace std;

class Point 
{
public:
	Point();
	Point(double x, double y, double z);
	~Point();
	double x, y, z;
	friend bool operator<(const Point & p_lhs, const Point & p_rhs);
	friend ostream& operator<<(ostream& os, const Point& point);
	friend bool operator==(const Point& lhs, const Point& rhs);
	friend bool operator!=(const Point& lhs, const Point& rhs);
	Point operator- (Point& left);
	Point operator+ (Point& left);
private:
};