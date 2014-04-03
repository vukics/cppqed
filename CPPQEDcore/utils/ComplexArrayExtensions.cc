// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ComplexArrayExtensions.tcc"

namespace blitzplusplus {

namespace dodirect {

using namespace linalg;
using namespace blitz::tensor;

template<> void doDirect<true >(CMatrix& m, const CVector& v1, const CVector& v2) {m=v1(i)*v2(j);}
template<> void doDirect<false>(CMatrix& m, const CVector& v1, const CVector& v2) {m=v1(i)+v2(j);}


} // dodirect


} // blitzplusplus
