// Copyright Raimar Sandner 2012–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-

#include "Core.h"

#include "LiouvilleanAveragedCommon.h"

#include <boost/python/exception_translator.hpp>


using namespace boost::python;

namespace pythonext {

namespace {
void pyInfiniteDetectedException(const structure::InfiniteDetectedException &)
{
  PyErr_SetString(PyExc_RuntimeError, "InfinateDetectedException in C++QED");
}
}

void export_io()
{
  register_exception_translator<structure::InfiniteDetectedException>(&pyInfiniteDetectedException);
}

} // pythonext