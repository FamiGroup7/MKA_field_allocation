#pragma once
#include "MkaUtils.h"
#include <iostream>

using namespace std;

class Point 
{
public:
	Point();
	Point(double x, double y, double z);
	Point(double x, double y, double z, int i);
	~Point();
	double x, y, z;
	int ind;
	friend bool operator<(const Point & p_lhs, const Point & p_rhs);
	friend ostream& operator<<(ostream& os, const Point& point);
	friend bool operator==(const Point& lhs, const Point& rhs);
	friend bool operator!=(const Point& lhs, const Point& rhs);
	Point operator- (Point& left);
	Point operator+ (Point& left);
	double & operator [] (size_t i);
	bool inside(double x0, double x1, double y0, double y1, double z0, double z1) const;
private:
};