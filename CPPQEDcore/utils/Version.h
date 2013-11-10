#ifndef UTILS_VERSION_H_INCLUEDED
#define UTILS_VERSION_H_INCLUEDED

#include <string>

extern std::string cppqed_versionstring;

std::string versionHelper();

void updateVersionstring(const std::string &s);

#endif // UTILS_VERSION_H_INCLUEDED