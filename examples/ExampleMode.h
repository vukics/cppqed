// -*- C++ -*-
/// [basic example mode]
#include "Free.h"
#include "TridiagonalHamiltonian.h"
#include "ElementLiouvillean.h"
#include "ElementAveraged.h"
/// [basic example mode]
#include "FreeExact.h"

/// [basic example mode]
using namespace structure;
using namespace free;


const Tridiagonal aop(size_t);
const Tridiagonal nop(size_t);


class PumpedLossyMode 
  : public Free, public TridiagonalHamiltonian<1,false>, public ElementLiouvillean<1,2>, public ElementAveraged<1>
{
public:
  typedef ElementAveraged<1>::LazyDensityOperator LazyDensityOperator;

  PumpedLossyMode(double delta, double kappa, dcomp eta, double n, size_t cutoff);

private:
  const Averages average_v(const LazyDensityOperator&) const;
  void           process_v(Averages&)                  const {}

};
/// [basic example mode]

class PumpedLossyModeIP
  : public Free, public FreeExact, public TridiagonalHamiltonian<1,true>, public ElementLiouvillean<1,2>, public ElementAveraged<1>
{
public:
  typedef ElementAveraged<1>::LazyDensityOperator LazyDensityOperator;

  PumpedLossyModeIP(double delta, double kappa, dcomp eta, double n, size_t cutoff);

  const dcomp get_z() const {return z_;}

private:
  void updateU(double dtDid) const;

  bool isUnitary_v() const {return true;}

  const Averages average_v(const LazyDensityOperator&) const;
  void           process_v(Averages&)                  const {}

  const dcomp z_; // Needed for updateU

};


const Tridiagonal aop(const PumpedLossyModeIP&);
