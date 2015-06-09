// Copyright András Vukics 2006–2015. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDELEMENTS_INTERACTIONS_MODECORRELATIONS_H_INCLUDED
#define   CPPQEDELEMENTS_INTERACTIONS_MODECORRELATIONS_H_INCLUDED

#include "ModeCorrelationsFwd.h"

#include "Mode_.h"

#include "ElementAveraged.h"
#include "Interaction.h"



class ModeCorrelations : public structure::ElementAveraged<2>
{
public:

  ModeCorrelations();

protected:
  typedef structure::ElementAveraged<2> EA_Base  ;
  typedef structure::       Averaged<1> Averaged1;

private:
  const Averages average_v(structure::NoTime, const LazyDensityOperator&) const;
  void           process_v(                   Averages&)                  const;

  const mode::AveragedQuadratures averagedMode_;

};



#endif // CPPQEDELEMENTS_INTERACTIONS_MODECORRELATIONS_H_INCLUDED
