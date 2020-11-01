// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Besides Evolution.h, BinarySystem gets included as well for scripts simulating binary quantum systems}
#ifndef CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTIONBINARY_H_INCLUDED
#define CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTIONBINARY_H_INCLUDED

#include "BinarySystem.h"

#include "Evolution.h"


template<typename SV_OR_DO>
inline
const auto
evolve(SV_OR_DO&& initial,
       binary::Ptr sys,
       const evolution::Pars<>& p,
       bool doStreaming=true, bool returnStreamedArray=false
      )
{
  return evolve<0>(std::forward<SV_OR_DO>(initial),sys,p,doStreaming,returnStreamedArray);
}


#endif // CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTIONBINARY_H_INCLUDED
