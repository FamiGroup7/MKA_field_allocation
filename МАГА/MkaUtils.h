#pragma once
#include <math.h>

const double EPS = 1e-15;

class MkaUtils
{
public:
	MkaUtils();
	~MkaUtils();
	static bool equals(double val1, double val2);
	static int compare(double val1, double val2);
};