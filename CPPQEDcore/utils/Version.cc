// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Version.h"

#include "core_version.h"

#include <blitz/bzconfig.h>
#include <gsl/gsl_version.h>
#include <boost/version.hpp>


std::string cppqed_versionstring = cppqed_core_version();

std::string versionHelper()
{
  return cppqed_versionstring + "\nAndras Vukics, vukics@users.sourceforge.net\n\nCompiled with\nBoost library collection : Version "+string(BOOST_LIB_VERSION).replace(1,1,".")+"\nGnu Scientific Library   : Version "+GSL_VERSION+"\nBlitz++ numerical library: "+BZ_PACKAGE_STRING+"\n\n";
}
