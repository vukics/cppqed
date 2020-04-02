// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDELEMENTS_INTERACTIONS_QBITMODECORRELATIONS_H_INCLUDED
#define   CPPQEDELEMENTS_INTERACTIONS_QBITMODECORRELATIONS_H_INCLUDED

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



#endif // CPPQEDELEMENTS_INTERACTIONS_QBITMODECORRELATIONS_H_INCLUDED
