// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_INTERACTIONS_MODECORRELATIONS_H_INCLUDED
#define   CPPQEDELEMENTS_INTERACTIONS_MODECORRELATIONS_H_INCLUDED

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
  const mode::Averages average_v(structure::NoTime, const quantumdata::LazyDensityOperator<2>&) const;
  void process_v(mode::Averages&) const;

  const mode::AveragedQuadratures averagedMode_;

};



#endif // CPPQEDELEMENTS_INTERACTIONS_MODECORRELATIONS_H_INCLUDED
