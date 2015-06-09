// Copyright András Vukics 2006–2015. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef UTILS_VERSION_H_INCLUEDED
#define UTILS_VERSION_H_INCLUEDED

#include <string>

extern std::string cppqed_versionstring;

std::string versionHelper();

void updateVersionstring(const std::string &s);

#endif // UTILS_VERSION_H_INCLUEDED