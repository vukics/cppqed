// -*- C++ -*-
#include "Free.h"
#include "TridiagonalHamiltonian.h"
#include "ElementLiouvillean.h"
#include "ElementAveraged.h"

using namespace structure;
using namespace free;


const Tridiagonal aop(size_t);
const Tridiagonal nop(size_t);


class PumpedLossyMode 
  : public Free, public TridiagonalHamiltonian<1,false>, public ElementLiouvillean<1,2>, public ElementAveraged<1>
{
public:
  typedef ElementAveraged<1>::LazyDensityOperator LazyDensityOperator;

  PumpedLossyMode(double omega, double kappa, dcomp eta, double n, size_t cutoff);

  const Averages average(const LazyDensityOperator&) const;
  void           process(Averages&)                  const {}

};
