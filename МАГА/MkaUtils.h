#pragma once
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <locale>
#include <bitset>
#include <Windows.h>
#include <time.h>
#include <chrono>
#include <sstream> // stringstream
#include <iomanip> // put_time
#include <string>  // string

const double EPS = 1e-15;

#define _USE_MATH_DEFINES

class MkaUtils
{
public:
	MkaUtils();
	~MkaUtils();
	static bool equals(double val1, double val2);
	static int compare(double val1, double val2);

	static std::string formattingTime(std::chrono::system_clock::time_point timestamp);
	static int Search_Binary(double arr[], int left, int right, double key);
};