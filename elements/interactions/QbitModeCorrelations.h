// -*- C++ -*-
#ifndef   ELEMENTS_INTERACTIONS_QBITMODECORRELATIONS_H_INCLUDED
#define   ELEMENTS_INTERACTIONS_QBITMODECORRELATIONS_H_INCLUDED

#include "ElementAveraged.h"


class QbitModeCorrelations : public structure::ElementAveraged<2>
{
public:

  QbitModeCorrelations();

protected:
  typedef structure::ElementAveraged<2> EA_Base  ;

private:
  const Averages average_v(const LazyDensityOperator&) const;
  void           process_v(Averages&)                  const {}

};



#endif // ELEMENTS_INTERACTIONS_QBITMODECORRELATIONS_H_INCLUDED
