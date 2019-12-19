// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDCORE_QUANTUMOPERATOR_SIGMAFWD_H_INCLUDED
#define   CPPQEDCORE_QUANTUMOPERATOR_SIGMAFWD_H_INCLUDED

namespace quantumoperator {


template<int L, int R>
// realizes the operator |L><R|
class Sigma;


template<int L, int R, typename OTHER, bool IS_HEAD>
class DirectProduct;


} // quantumoperator

#endif // CPPQEDCORE_QUANTUMOPERATOR_SIGMAFWD_H_INCLUDED
