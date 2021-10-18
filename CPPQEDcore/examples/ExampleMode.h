// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// [basic example mode]
#include "Free.h"
#include "ElementLiouvillean.h"
#include "ElementAveraged.h"

#include "TridiagonalHamiltonian.h"

using quantumoperator::TridiagonalHamiltonian, structure::NoTime, structure::OneTime;
using namespace structure::freesystem;


const Tridiagonal aop(size_t dim); // ladder operator of the given dimension
const Tridiagonal nop(size_t dim); // number operator "


namespace basic {


class PumpedLossyMode 
  : public structure::Free, public TridiagonalHamiltonian<1,false>, public structure::ElementLiouvilleanStrategies<1,2>, public structure::ElementAveraged<1>
{
public:
  PumpedLossyMode(double delta, double kappa, dcomp eta, double nTh, size_t cutoff);

private:
  const Averages average_v(NoTime, const structure::freesystem::LazyDensityOperator&) const;

};


} // basic

inline const Tridiagonal aop(const basic::PumpedLossyMode& mode) {return aop(mode.getDimension());} // just for convenience

/// [basic example mode]

#include "FreeExact.h" // Further inculdes and using directives the same as above

namespace basic {


class PumpedLossyModeIP
  : public structure::Free, public structure::FreeExact<false>, public TridiagonalHamiltonian<1,true>,
    public structure::ElementLiouvilleanStrategies<1,2>, public structure::ElementAveraged<1>
{
public:
  PumpedLossyModeIP(double delta, double kappa, dcomp eta, double nTh, size_t cutoff);

  const dcomp get_z() const {return z_;}

private:
  void updateU(OneTime) const override;

  bool applicableInMaster_v() const override {return true;}

  const Averages average_v(NoTime, const structure::freesystem::LazyDensityOperator&) const override;

  const dcomp z_; // Needed for updateU

};


} // basic


namespace hierarchical {

// All inculdes and using directives the same as above

class ModeBase : public structure::Free, public structure::ElementLiouvillean<1,2>, public structure::ElementAveraged<1>
{
public:
  typedef std::shared_ptr<ModeBase> Ptr;
  
  using StateVectorLow = ::quantumdata::StateVectorLow<1>;
  
  using LazyDensityOperator = ::quantumdata::LazyDensityOperator<1>;

protected:
  ModeBase(double kappa, double nTh, size_t cutoff);
  
private:
  void doActWithJ(NoTime, StateVectorLow&, LindbladNo<0>) const override;
  void doActWithJ(NoTime, StateVectorLow&, LindbladNo<1>) const override;
  
  double rate(NoTime, const LazyDensityOperator&, LindbladNo<0>) const override;
  double rate(NoTime, const LazyDensityOperator&, LindbladNo<1>) const override;
  
  const Averages average_v(NoTime, const LazyDensityOperator&) const override;

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
  : public ModeBase, public structure::FreeExact<false>, public TridiagonalHamiltonian<1,true >
{
public:
  PumpedLossyModeIP(double delta, double kappa, dcomp eta, double nTh, size_t cutoff);

  const dcomp get_z() const {return z_;}

private:
  void updateU(OneTime) const override;

  bool applicableInMaster_v() const override {return true;}

  const dcomp z_; // Needed for updateU

};


} // hierarchical


const Tridiagonal aop(const hierarchical::ModeBase& mode);
