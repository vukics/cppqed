// -*- C++ -*-
/// [basic example mode]
#include "Free.h"
#include "ElementLiouvillean.h"
#include "ElementAveraged.h"

#include "TridiagonalHamiltonian.h"

using namespace structure;
using quantumoperator::TridiagonalHamiltonian;
using namespace freesystem;


const Tridiagonal aop(size_t dim); // ladder operator of the given dimension
const Tridiagonal nop(size_t dim); // number operator "


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


} // basic

inline const Tridiagonal aop(const basic::PumpedLossyMode& mode) {return aop(mode.getDimension());} // just for convenience

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

// All inculdes and using directives the same as above

class ModeBase : public Free, public ElementLiouvillean<1,2>, public ElementAveraged<1>
{
public:
  typedef ElementLiouvillean<1,2>::StateVectorLow StateVectorLow     ;
  typedef ElementAveraged<1>::LazyDensityOperator LazyDensityOperator;

protected:
  ModeBase(double kappa, double nTh, size_t cutoff);
  
private:
  void doActWithJ(NoTime, StateVectorLow&, LindbladNo<0>) const;
  void doActWithJ(NoTime, StateVectorLow&, LindbladNo<1>) const;
  
  double rate(NoTime, const LazyDensityOperator&, LindbladNo<0>) const;
  double rate(NoTime, const LazyDensityOperator&, LindbladNo<1>) const;
  
  const Averages average_v(NoTime, const LazyDensityOperator&) const;

  const double kappa_, nTh_; // needed for calculating jumps & rates 
  
};

} // hierarchical


namespace hierarchical {


class PumpedLossyMode
  : public ModeBase, public TridiagonalHamiltonian<1,false>
{
public:
  PumpedLossyMode(double delta, double kappa, dcomp eta, double nTh, size_t cutoff);

};


class PumpedLossyModeIP
  : public ModeBase, public FreeExact<false>, public TridiagonalHamiltonian<1,true >
{
public:
  PumpedLossyModeIP(double delta, double kappa, dcomp eta, double nTh, size_t cutoff);

  const dcomp get_z() const {return z_;}

private:
  void updateU(OneTime) const;

  bool isUnitary_v() const {return true;}

  const dcomp z_; // Needed for updateU

};


} // hierarchical


const Tridiagonal aop(const hierarchical::ModeBase& mode);
