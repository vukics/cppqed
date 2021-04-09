// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ComplexArrayExtensions.h"

using namespace linalg;
using namespace blitz::tensor;

template<> void blitzplusplus::dodirect::_<true >(CMatrix& m, const CVector& v1, const CVector& v2) {m=v1(i)*v2(j);}
template<> void blitzplusplus::dodirect::_<false>(CMatrix& m, const CVector& v1, const CVector& v2) {m=v1(i)+v2(j);}
