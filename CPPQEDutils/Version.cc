// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Version.h"

#include <blitz/config.h>
#include <gsl/gsl_version.h>
#include <boost/version.hpp>

using namespace std;

string versionHelper()
{
  return "http://github.com/vukics/cppqed\ncommit# " +
         string(cppqed_GIT_SHA1) + 
         "\n\nCompiled with\nBoost library collection : Version " +
         string(BOOST_LIB_VERSION).replace(1,1,".") +
         "\nGnu Scientific Library   : Version " +
         GSL_VERSION +
         "\nBlitz++ numerical library: http://github.com/vukics/blitz #" + string(blitz_GIT_SHA1)+"\n\n";
}
