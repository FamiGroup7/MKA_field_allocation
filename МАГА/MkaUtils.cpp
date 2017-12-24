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

std::string MkaUtils::formattingTime(std::chrono::system_clock::time_point timestamp)
{
	auto hours = std::chrono::time_point_cast<std::chrono::hours>(timestamp);
	auto minutes = std::chrono::time_point_cast<std::chrono::minutes>(timestamp);
	auto seconds = std::chrono::time_point_cast<std::chrono::seconds>(timestamp);

	auto frminutes = std::chrono::duration_cast<std::chrono::minutes>(timestamp - hours);
	auto frsec = std::chrono::duration_cast<std::chrono::seconds>(timestamp - minutes);
	auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(timestamp - seconds);

	char buff[100];
	snprintf(buff, sizeof(buff), "%02d:%02d:%02dM%03d", frminutes.count(), frsec.count(), milliseconds.count());
	return std::string(buff);
}
