// -*- C++ -*-
#ifndef   SIGMA_TYPE_OPERATOR_FWD_INCLUDED
#define   SIGMA_TYPE_OPERATOR_FWD_INCLUDED

namespace quantumoperator {


template<int L, int R>
// realizes the operator |L><R|
class Sigma;


template<int L, int R, typename OTHER, bool IS_HEAD>
class DirectProduct;


} // quantumoperator

#endif // SIGMA_TYPE_OPERATOR_FWD_INCLUDED
