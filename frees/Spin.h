// -*- C++ -*-
#ifndef   ELEMENTS_FREES_SPIN_H_INCLUDED
#define   ELEMENTS_FREES_SPIN_H_INCLUDED

#include "SpinFwd.h"

#include "ElementLiouvillean.h"
#include "ElementAveraged.h"
#include "Free.h"
#include "FreeExact.h"
#include "TridiagonalHamiltonian.h"

#include "ParsFwd.h"
#include "SmartPtr.h"


// A general Spin yet incomplete
// Note: jump is not yet implemented, only "Hamiltonian" decay.

namespace spin {

using namespace structure::free;


typedef boost::shared_ptr<const SpinBase> Ptr;


const Tridiagonal splus(Ptr);

inline const Tridiagonal sminus(Ptr spin) {return splus(spin).dagger();}

inline const Tridiagonal sx(Ptr spin) {return (splus(spin)+sminus(spin))/2;}
inline const Tridiagonal sy(Ptr spin) {return (splus(spin)-sminus(spin))/(2.*DCOMP_I);}

const Tridiagonal sn(Ptr);


const Tridiagonal sz(Ptr);


struct Pars
{
  size_t &twoS, &dim;
  double &theta, &phi;
  double &omega, &gamma;

  Pars(parameters::ParameterTable&, const std::string& ="");

};


} // spin


class SpinBase
  : public structure::Free, 
    public structure::ElementAveraged<1>
{
public:
  typedef structure::ElementAveraged<1>::LazyDensityOperator LazyDensityOperator;

  SpinBase(size_t twoS, double theta, double phi, double omega, double gamma, size_t dim=0);

  size_t getTwoS() const {return twoS_;}

  double getTheta() const {return theta_;}
  double getPhi  () const {return   phi_;}

  double getOmega() const {return omega_;}
  double getGamma() const {return gamma_;}

  const dcomp get_z() const {return dcomp(gamma_,omega_);}
  // This plays analogous role as Z in Mode

private:
  void process_v(Averages&) const;

  const Averages average_v(const LazyDensityOperator&) const;

  const size_t twoS_;

  const double theta_, phi_, omega_, gamma_, s_;

};

// NOTE if p.gamma is nonzero, Spin & SpinSch will contain the „Hamiltonian” part of the loss!

class Spin 
  : public SpinBase,
    public structure::FreeExact<false>
{
public:

  Spin(const spin::Pars& p) 
    : SpinBase(p.twoS,0.,0.,p.omega,p.gamma,p.dim), FreeExact(getTotalDimension())
      // here, we need to assume that the Hamiltonian is diagonal in the given basis, that is why we set theta=phi=0.
  {}

private:
  void updateU(double) const;

  bool isUnitary_v() const {return !getGamma();}

};


class SpinSch
  : public SpinBase,
    public structure::TridiagonalHamiltonian<1,false>
{
public:

  SpinSch(const spin::Pars& p) 
    : SpinBase(p.twoS,p.theta,p.phi,p.omega,p.gamma,p.dim),
      structure::TridiagonalHamiltonian<1,false>(-get_z()*spin::sn(cpputils::nonOwningConstSharedPtr(this)))
  {
    getParsStream()<<"# Schrodinger picture."<<std::endl;
  }

};


class LossySpin
  : public Spin,
    public structure::ElementLiouvillean<1,1,false>
{
public:
  
  LossySpin(const spin::Pars& p) : Spin(p), structure::ElementLiouvillean<1,1,false>("Spin","S-") {}
  
private:
  
  void   doActWithJ (      structure::free::StateVectorLow     &) const;
  double probability(const structure::free::LazyDensityOperator&) const {return -1.;}
  
};


#endif // ELEMENTS_FREES_SPIN_H_INCLUDED
