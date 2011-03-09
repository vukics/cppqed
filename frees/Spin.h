// -*- C++ -*-
#ifndef   _SPIN_ELEMENT_INCLUDED
#define   _SPIN_ELEMENT_INCLUDED

#include "SpinFwd.h"

#include "ElementLiouvillean.h"
#include "ElementAveraged.h"
#include "Free.h"
#include "FreeExact.h"
#include "TridiagonalHamiltonian.h"

#include "ParsFwd.h"

// A general Spin yet incomplete
// Note: jump is not yet implemented, only "Hamiltonian" decay.

namespace spin {

using namespace structure::free;


const Tridiagonal splus (size_t twos, size_t dim=0);

inline const Tridiagonal sminus(size_t twos, size_t dim=0) {return splus(twos,dim).dagger();}

inline const Tridiagonal sx(size_t twos, size_t dim=0) {return (splus(twos,dim)+sminus(twos,dim))/2;}
inline const Tridiagonal sy(size_t twos, size_t dim=0) {return (splus(twos,dim)-sminus(twos,dim))/(2.*DCOMP_I);}

const Tridiagonal sz(size_t, size_t dim=0);

const Frequencies freqs(const SpinBase*);


struct Pars
{
  size_t &twos, &dim;
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

  SpinBase(size_t twos, double omega, double gamma, size_t dim=0);

  double getOmega() const {return omega_;}
  double getGamma() const {return gamma_;}

  const dcomp get_z() const {return dcomp(gamma_,omega_);}
  // This plays analogous role as Z in Mode

private:
  void process(Averages&) const;

  const Averages average(const LazyDensityOperator&) const;

  double omega_, gamma_;

};



class Spin 
  : public SpinBase,
    public structure::FreeExact
{
public:

  Spin(const spin::Pars& p) 
    : SpinBase(p.twos,p.omega,p.gamma,p.dim), FreeExact(getTotalDimension())
  {}

private:
  bool isUnitary() const {return !getGamma();}

  void updateU(double) const;

};


class SpinSch
  : public SpinBase,
    public structure::TridiagonalHamiltonian<1,false>
{
public:

  SpinSch(const spin::Pars& p) 
    : SpinBase(p.twos,p.omega,p.gamma,p.dim),
      structure::TridiagonalHamiltonian<1,false>(-get_z()*spin::sz(p.twos))
  {
    getParsStream()<<"# Schrodinger picture."<<std::endl;
  }

};


#endif // _SPIN_ELEMENT_INCLUDED
