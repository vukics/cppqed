// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines the spin-bundle}
#ifndef   CPPQEDELEMENTS_FREES_SPIN_H_INCLUDED
#define   CPPQEDELEMENTS_FREES_SPIN_H_INCLUDED

#include "ElementLiouvillean.h"
#include "ElementAveraged.h"
#include "Free.h"
#include "FreeExact.h"
#include "TridiagonalHamiltonian.h"

#include "Pars.h"


class SpinBase;

// A general Spin yet incomplete
// Implements the following Hamiltonian: 
// (omega-i*gamma)*S_z
// with the usual Liouvillean

namespace spin {

using namespace structure::freesystem;


typedef std::shared_ptr<const SpinBase> Ptr;

Tridiagonal splus(size_t dim, size_t twoS, dcomp z, bool isExact);

inline auto sminus(size_t dim, size_t twoS, dcomp z, bool isExact) {return splus(dim,twoS,z,isExact).dagger();}
inline auto sx(size_t dim, size_t twoS, dcomp z, bool isExact) {return (splus(dim,twoS,z,isExact)+sminus(dim,twoS,z,isExact))/2;}
inline auto sy(size_t dim, size_t twoS, dcomp z, bool isExact) {return (splus(dim,twoS,z,isExact)-sminus(dim,twoS,z,isExact))/(2.*DCOMP_I);}


Tridiagonal splus(Ptr spin);

inline auto sminus(Ptr spin) {return splus(spin).dagger();}
inline auto sx(Ptr spin) {return (splus(spin)+sminus(spin))/2;}
inline auto sy(Ptr spin) {return (splus(spin)-sminus(spin))/(2.*DCOMP_I);}

Tridiagonal sz(size_t dim, size_t twoS);

Tridiagonal sz(Ptr spin);

inline auto sn(size_t dim, size_t twoS, dcomp z, bool isExact, double theta, double phi)
{
  return sin(theta)*(cos(phi)*sx(dim,twoS,z,isExact)+sin(phi)*sy(dim,twoS,z,isExact))+cos(theta)*sz(dim,twoS);
}

Tridiagonal sn(Ptr);


struct Pars
{
  size_t &twoS, &dim;
  double &theta, &phi;
  double &omega, &gamma;

  Pars(parameters::Table&, const std::string& ="");

};


class Liouvillean
  : public structure::ElementLiouvillean<1,1,false>
{
protected:
  Liouvillean(size_t twoS, double gamma) : structure::ElementLiouvillean<1,1,false>("Spin","S-"), twoS_(twoS), gamma_(gamma) {}
  
private:
  void   doActWithJ(structure::NoTime,       structure::freesystem::StateVectorLow     &) const ;
  double rate      (structure::NoTime, const structure::freesystem::LazyDensityOperator& matrix) const ; // {return -1.;}

  const size_t twoS_;
  const double gamma_;
  
};


} // spin


/** \todo Implement some spherical-coordinates class and use it here */
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

  const dcomp get_z() const {return dcomp(gamma_,omega_);} ///< This plays analogous role to \f$z\f$ in Mode

private:
  void process_v(spin::Averages&) const;

  const spin::Averages average_v(structure::NoTime, const LazyDensityOperator&) const;

  const size_t twoS_;

  const double theta_, phi_, omega_, gamma_, s_;

};

// NOTE if p.gamma is nonzero, Spin & SpinSch will contain the “Hamiltonian” part of the loss!

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
  void updateU(structure::OneTime) const;

  bool applicableInMaster_v() const {return !getGamma();}

};


class SpinSch
  : public SpinBase,
    public quantumoperator::TridiagonalHamiltonian<1,false>
{
public:

  SpinSch(const spin::Pars& p) 
    : SpinBase(p.twoS,p.theta,p.phi,p.omega,p.gamma,p.dim),
      quantumoperator::TridiagonalHamiltonian<1,false>(-get_z()*spin::sn(getDimension(),p.twoS,get_z(),false,p.theta,p.phi))
  {
    getParsStream()<<"Schrodinger picture."<<std::endl;
  }

};


class LossySpin
  : public Spin,
    public spin::Liouvillean
{
public:
  
  LossySpin(const spin::Pars& p) : Spin(p), spin::Liouvillean(p.twoS,p.gamma) {}
  
};


class LossySpinSch
  : public SpinSch,
    public spin::Liouvillean
{
public:
  
  LossySpinSch(const spin::Pars& p) : SpinSch(p), spin::Liouvillean(p.twoS,p.gamma) {}
  
};


#endif // CPPQEDELEMENTS_FREES_SPIN_H_INCLUDED
