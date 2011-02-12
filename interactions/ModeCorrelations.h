// -*- C++ -*-
#ifndef   MODE_CORRELATIONS_INCLUDED
#define   MODE_CORRELATIONS_INCLUDED

#include "ModeCorrelationsFwd.h"

#include "Mode.h"

#include "ElementAveraged.h"
#include "Interaction.h"



class ModeCorrelations : public structure::ElementAveraged<2>
{
protected:
  typedef structure::ElementAveraged<2> EA_Base;

  ModeCorrelations();

private:
  const Averages average(const LazyDensityOperator&) const;// {return Averages(18);}
  void           process(Averages&)                  const;// {}

  const mode::AveragedQuadratures averagedMode_;

};



#endif // MODE_CORRELATIONS_INCLUDED
