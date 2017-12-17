#pragma once
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <locale>
#include <bitset>
#include <Windows.h>

const double EPS = 1e-15;

class MkaUtils
{
public:
	MkaUtils();
	~MkaUtils();
	static bool equals(double val1, double val2);
	static int compare(double val1, double val2);
};