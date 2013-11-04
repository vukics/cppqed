// -*- C++ -*-
#ifndef   INTERACTIONS_MODECORRELATIONS_H_INCLUDED
#define   INTERACTIONS_MODECORRELATIONS_H_INCLUDED

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



#endif // INTERACTIONS_MODECORRELATIONS_H_INCLUDED