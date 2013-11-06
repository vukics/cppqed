// -*- C++ -*-
#ifndef   INTERACTIONS_QBITMODECORRELATIONS_H_INCLUDED
#define   INTERACTIONS_QBITMODECORRELATIONS_H_INCLUDED

#include "ElementAveraged.h"


class QbitModeCorrelations : public structure::ElementAveraged<2>
{
public:

  QbitModeCorrelations();

protected:
  typedef structure::ElementAveraged<2> EA_Base;
  typedef structure::ElementAveraged<2>::Time NoTime;

private:
  const Averages average_v(NoTime, const LazyDensityOperator&) const;

};



#endif // INTERACTIONS_QBITMODECORRELATIONS_H_INCLUDED
