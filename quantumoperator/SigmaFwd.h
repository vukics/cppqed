// -*- C++ -*-
#ifndef   QUANTUMOPERATOR_SIGMAFWD_H_INCLUDED
#define   QUANTUMOPERATOR_SIGMAFWD_H_INCLUDED

namespace quantumoperator {


template<int L, int R>
// realizes the operator |L><R|
class Sigma;


template<int L, int R, typename OTHER, bool IS_HEAD>
class DirectProduct;


} // quantumoperator

#endif // QUANTUMOPERATOR_SIGMAFWD_H_INCLUDED
