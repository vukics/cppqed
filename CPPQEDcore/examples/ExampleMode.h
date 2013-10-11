// -*- C++ -*-
/// [basic example mode]
#include "Free.h"
#include "TridiagonalHamiltonian.h"
#include "ElementLiouvillean.h"
#include "ElementAveraged.h"

using namespace structure;
using namespace free;


const Tridiagonal aop(size_t);
const Tridiagonal nop(size_t);


namespace basic {


class PumpedLossyMode 
  : public Free, public TridiagonalHamiltonian<1,false>, public ElementLiouvilleanStrategies<1,2>, public ElementAveraged<1>
{
public:
  typedef ElementAveraged<1>::LazyDensityOperator LazyDensityOperator;

  PumpedLossyMode(double delta, double kappa, dcomp eta, double nTh, size_t cutoff);

private:
  const Averages average_v(NoTime, const LazyDensityOperator&) const;

};

} //basic

/// [basic example mode]

#include "FreeExact.h" // Further inculdes and using directives the same as above

namespace basic {


class PumpedLossyModeIP
  : public Free, public FreeExact<false>, public TridiagonalHamiltonian<1,true>, public ElementLiouvilleanStrategies<1,2>, public ElementAveraged<1>
{
public:
  typedef ElementAveraged<1>::LazyDensityOperator LazyDensityOperator;

  PumpedLossyModeIP(double delta, double kappa, dcomp eta, double nTh, size_t cutoff);

  const dcomp get_z() const {return z_;}

private:
  void updateU(OneTime) const;

  bool isUnitary_v() const {return true;}

  const Averages average_v(NoTime, const LazyDensityOperator&) const;

  const dcomp z_; // Needed for updateU

};


} // basic



namespace hierarchical {


class ModeBase : public Free, public ElementLiouvillean<1,2>, public ElementAveraged<1>
{
public:
  typedef ElementLiouvillean<1,2>::StateVectorLow StateVectorLow     ;
  typedef ElementAveraged<1>::LazyDensityOperator LazyDensityOperator;

protected:
  ModeBase(double kappa, double nTh, size_t cutoff);
  
private:
  void doActWithJ(NoTime, StateVectorLow&, JumpNo<0>) const;
  void doActWithJ(NoTime, StateVectorLow&, JumpNo<1>) const;
  
  double rate(NoTime, const LazyDensityOperator&, JumpNo<0>) const;
  double rate(NoTime, const LazyDensityOperator&, JumpNo<1>) const;
  
  const Averages average_v(NoTime, const LazyDensityOperator&) const;

  const double kappa_, nTh_;
  
};


class PumpedLossyMode
  : public ModeBase, public TridiagonalHamiltonian<1,false>
{
public:
  PumpedLossyMode(double delta, double kappa, dcomp eta, double n, size_t cutoff);

};


class PumpedLossyModeIP
  : public ModeBase, public TridiagonalHamiltonian<1,true >
{
public:
  PumpedLossyModeIP(double delta, double kappa, dcomp eta, double n, size_t cutoff);

  const dcomp get_z() const {return z_;}

private:
  void updateU(OneTime) const;

  bool isUnitary_v() const {return true;}

  const dcomp z_; // Needed for updateU

};


const Tridiagonal aop(const ModeBase& mode);

} // hierarchical