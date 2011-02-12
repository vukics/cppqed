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


const Tridiagonal splus (size_t twos);

inline const Tridiagonal sminus(size_t twos) {return splus(twos).dagger();}

inline const Tridiagonal sx(size_t twos) {return (splus(twos)+sminus(twos))/2;}
inline const Tridiagonal sy(size_t twos) {return (splus(twos)-sminus(twos))/(2.*DCOMP_I);}

const Tridiagonal sz(size_t);

const Frequencies freqs(const SpinBase*);


struct Pars
{
  size_t &twos;
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

  SpinBase(size_t twos, double omega, double gamma);

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
    : SpinBase(p.twos,p.omega,p.gamma), FreeExact(getTotalDimension())
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
    : SpinBase(p.twos,p.omega,p.gamma),
      structure::TridiagonalHamiltonian<1,false>(-get_z()*spin::sz(p.twos))
  {}

};


#endif // _SPIN_ELEMENT_INCLUDED
