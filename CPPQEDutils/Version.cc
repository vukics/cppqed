// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Version.h"

#include <boost/version.hpp>

#include <Eigen/Eigen>

using namespace std;

string versionHelper()
{
  return "http://github.com/vukics/cppqed\ncommit# " +
         string(cppqed_GIT_SHA1) +
#ifndef NDEBUG
         "\nDEBUG build" +
#endif // NDEBUG
         "\n\nCompiled with\nBoost library collection : Version " +
         string(BOOST_LIB_VERSION).replace(1,1,".") +
         "\nEigen library            : Version "+to_string(EIGEN_WORLD_VERSION)+"."+to_string(EIGEN_MAJOR_VERSION)+"."+to_string(EIGEN_MINOR_VERSION)+
         "\n\n";         
        
}
