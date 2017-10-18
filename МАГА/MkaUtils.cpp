#include "MkaUtils.h"

MkaUtils::MkaUtils()
{
}

MkaUtils::~MkaUtils()
{
}

bool MkaUtils::equals(double val1, double val2) {
	if (fabs(val1 - val2) < EPS)return true;
	return false;
}

int MkaUtils::compare(double val1, double val2) {
	if (equals(val1, val2)) return 0;
	double diff = val1 - val2;
	return (diff > 0) - (diff < 0);
}