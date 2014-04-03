// Copyright Raimar Sandner 2012–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef CUSTOMELEMENTSEXAMPLE_INTERACTIONS_MOMENTUMCORRELATION_H_INCLUDED
#define CUSTOMELEMENTSEXAMPLE_INTERACTIONS_MOMENTUMCORRELATION_H_INCLUDED

#include "ElementAveraged.h"
#include "Interaction.h"
#include "Particle.h"
#include "LazyDensityOperator.h"

class MomentumCorrelation : public structure::Interaction<2>, public structure::ElementAveraged<2>
{
public:
  typedef structure::ElementAveraged<2> EA_Base;
  typedef quantumdata::LazyDensityOperator<2> LazyDensityOperator;
  MomentumCorrelation(particle::Ptr, particle::Ptr);
  
private:
  const Averages average_v(Time t, const LazyDensityOperator&) const;
  void           process_v(Averages&)                  const {}
};

#endif // CUSTOMELEMENTSEXAMPLE_INTERACTIONS_MOMENTUMCORRELATION_H_INCLUDED
