#include "ComplexArrayExtensions.tcc"

namespace blitzplusplus {

namespace dodirect {

using namespace linalg;
using namespace blitz::tensor;

template<> void doDirect<true >(CMatrix& m, const CVector& v1, const CVector& v2) {m=v1(i)*v2(j);}
template<> void doDirect<false>(CMatrix& m, const CVector& v1, const CVector& v2) {m=v1(i)+v2(j);}


} // dodirect


} // blitzplusplus
