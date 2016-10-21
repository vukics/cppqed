// Copyright András Vukics 2006–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDELEMENTS_FREES_MULTILEVEL_FWD_H_INCLUDED
#define   CPPQEDELEMENTS_FREES_MULTILEVEL_FWD_H_INCLUDED

#include "AveragingUtilsFwd.h"

/*

Convention is the following: an (l,m) pair in VP (or VL) below means a sigma_lm=|l><m|

This yields an 

H_elementaryPumping=i*(conj(eta)*sigma-eta*sigmadagger)

term in the Hamiltonian

If a pair (l,m) is in VL, this means a jump with jump operator
sigma_lm. This means that the rate of the jump is proportional
to 

<sigmadagger_lm*sigma_lm>=<|m><l|*|l><m|>=<|m><m|>=rho_mm

*/

namespace multilevel {

struct ElementaryComplexFreqs;


} // multilevel


template<int NL> 
// NL is the number of levels
class MultiLevelBase;


template<int NL, typename Averaged=ReducedDensityOperator<1> >
class MultiLevel;

template<int NL, typename Averaged=ReducedDensityOperator<1> >
class MultiLevelSch;


template<int NL, typename VP, typename Averaged=ReducedDensityOperator<1> > 
// VP is a compile-time container of pairs, specifying which
// transitions are pumped It should model a Boost.Fusion sequence,
// which stores the pairs for compile-time use and the pump Rabi
// frequencies for run-time use.
class PumpedMultiLevel;

template<int NL, typename VP, typename Averaged=ReducedDensityOperator<1> >
class PumpedMultiLevelSch;


template<int NL, typename VL, typename Averaged=ReducedDensityOperator<1> > 
// VL is a compile-time container of pairs, specifying which transitions have radiative loss
class LossyMultiLevel;

template<int NL, typename VL, typename Averaged=ReducedDensityOperator<1> >
class LossyMultiLevelSch;

template<int NL, typename VL, typename Averaged=ReducedDensityOperator<1> >
class LossyMultiLevelUIP;


template<int NL, typename VP, typename VL, typename Averaged=ReducedDensityOperator<1> >
class PumpedLossyMultiLevel;

template<int NL, typename VP, typename VL, typename Averaged=ReducedDensityOperator<1> >
class PumpedLossyMultiLevelSch;

template<int NL, typename VP, typename VL, typename Averaged=ReducedDensityOperator<1> >
class PumpedLossyMultiLevelUIP;


namespace multilevel {

template<int,int>
class Pump;

template<int,int>
class Decay;


struct Pars;
struct ParsPumped;
struct ParsLossy;

template<int NL, typename VP, typename VL>
struct ParsPumpedLossy;

// struct PrepError;

} // multilevel


#endif // CPPQEDELEMENTS_FREES_MULTILEVEL_FWD_H_INCLUDED
