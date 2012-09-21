// -*- C++ -*-
#ifndef   ELEMENTS_INTERACTIONS_MODECORRELATIONS_H_INCLUDED
#define   ELEMENTS_INTERACTIONS_MODECORRELATIONS_H_INCLUDED

#include "ModeCorrelationsFwd.h"

#include "Mode_.h"

#include "ElementAveraged.h"
#include "Interaction.h"



class ModeCorrelations : public structure::ElementAveraged<2>
{
protected:
  typedef structure::ElementAveraged<2> EA_Base  ;
  typedef structure::       Averaged<1> Averaged1;

  ModeCorrelations();

private:
  const Averages average(const LazyDensityOperator&) const;
  void           process(Averages&)                  const;

  const mode::AveragedQuadratures averagedMode_;

};



#endif // ELEMENTS_INTERACTIONS_MODECORRELATIONS_H_INCLUDED
