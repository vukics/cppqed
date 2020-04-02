// Copyright Raimar Sandner 2012–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-

#ifndef CPYPYQED_CORE_CORE_H_INCLUDED
#define CPYPYQED_CORE_CORE_H_INCLUDED

#include "PythonExtension.h"
#include "Namespaces.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_ARRAY_UNIQUE_SYMBOL core_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>

#endif // CPYPYQED_CORE_CORE_H_INCLUDED
