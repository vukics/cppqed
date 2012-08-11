// -*- C++ -*-
#ifndef   MULTI_LEVEL_SYSTEM_FWD_INCLUDED
#define   MULTI_LEVEL_SYSTEM_FWD_INCLUDED

#include "ElementAveragedFwd.h"

/*

Convention is the following: an (l,m) pair in VP (or VL) below means a sigma_lm=|l><m|

This yields an 

H_elementaryPumping=i*(conj(eta)*sigma-eta*sigmadagger)

term in the Hamiltonian

If a pair (l,m) is in VL, this means a jump with jump operator
sigma_lm. This means that the probability of the jump is proportional
to 

<sigmadagger_lm*sigma_lm>=<|m><l|*|l><m|>=<|m><m|>=rho_mm

*/

namespace multilevel {

using structure::averaged::DiagonalDO;

struct ElementaryComplexFreqs;


} // multilevel


template<int NL> 
// NL is the number of levels
class MultiLevelBase;


template<int NL, typename Averaged=multilevel::DiagonalDO>
class MultiLevel;

template<int NL, typename Averaged=multilevel::DiagonalDO>
class MultiLevelSch;


template<int NL, typename VP, typename Averaged=multilevel::DiagonalDO> 
// VP is a compile-time container of pairs, specifying which
// transitions are pumped It should model a Boost.Fusion sequence,
// which stores the pairs for compile-time use and the pump Rabi
// frequencies for run-time use.
class PumpedMultiLevel;

template<int NL, typename VP, typename Averaged=multilevel::DiagonalDO>
class PumpedMultiLevelSch;


template<int NL, typename VL, typename Averaged=multilevel::DiagonalDO> 
// VL is a compile-time container of pairs, specifying which transitions have radiative loss
class LossyMultiLevel;

template<int NL, typename VL, typename Averaged=multilevel::DiagonalDO>
class LossyMultiLevelSch;

template<int NL, typename VL, typename Averaged=multilevel::DiagonalDO>
class LossyMultiLevelUIP;


template<int NL, typename VP, typename VL, typename Averaged=multilevel::DiagonalDO>
class PumpedLossyMultiLevel;

template<int NL, typename VP, typename VL, typename Averaged=multilevel::DiagonalDO>
class PumpedLossyMultiLevelSch;

template<int NL, typename VP, typename VL, typename Averaged=multilevel::DiagonalDO>
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


#endif // MULTI_LEVEL_SYSTEM_FWD_INCLUDED
