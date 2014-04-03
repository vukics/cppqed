// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Version.h"

#include "core_version.h"

#include <blitz/bzconfig.h>
#include <gsl/gsl_version.h>
#include <boost/version.hpp>


std::string cppqed_versionstring = "# " + cppqed_core_version();

std::string versionHelper()
{
  return cppqed_versionstring + "\n# Andras Vukics, vukics@users.sourceforge.net\n\n# Compiled with\n# Boost library collection : Version "+string(BOOST_LIB_VERSION).replace(1,1,".")+"\n# Gnu Scientific Library   : Version "+GSL_VERSION+"\n# Blitz++ numerical library: Config date: "+BZ__config_date+"\n\n";
}

void updateVersionstring(const std::string &s)
{
  cppqed_versionstring=s;
}