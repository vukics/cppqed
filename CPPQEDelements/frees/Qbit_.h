// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines the qbit-bundle (tackling the dynamics of a single qbit)}
#ifndef CPPQEDELEMENTS_FREES_QBIT__H_INCLUDED
#define CPPQEDELEMENTS_FREES_QBIT__H_INCLUDED

#include "ParsQbit.h"

#include "Mode_.h"
#include "StateVector.h"

// Qbit traced back to Mode with dim=2

namespace qbit {

using namespace structure::freesystem;

using mode::NoTime; using mode::Averages;

const std::string keyTitle="Qbit";


class Averaged
  : public structure::ElementAveraged<1>
{
public:
  typedef structure::ElementAveraged<1> Base;

  Averaged();

private:
  const Averages average_v(NoTime, const qbit::LazyDensityOperator&) const;

};


} // qbit




class QbitBase 
  : public ModeBase, public qbit::Averaged
{
public:
  explicit QbitBase(const RealFreqs& ={}, const ComplexFreqs& ={});

  explicit QbitBase(const ComplexFreqs& complexFreqs) : QbitBase({},complexFreqs) {}
  explicit QbitBase(RealFreqsInitializer rf, ComplexFreqsInitializer cf={}) : QbitBase(RealFreqs(rf),ComplexFreqs(cf)) {}
  explicit QbitBase(ComplexFreqsInitializer cf) : QbitBase({},cf) {}
  explicit QbitBase(RF rf, CF cf=CF()) : QbitBase(RealFreqsInitializer{rf}, cf==CF() ? ComplexFreqsInitializer{} : ComplexFreqsInitializer{cf}) {}
  explicit QbitBase(CF cf) : QbitBase(ComplexFreqsInitializer{cf}) {}
  explicit QbitBase(RealFreqsInitializer rf, CF cf) : QbitBase(rf,{cf}) {}
  explicit QbitBase(RF rf, ComplexFreqsInitializer cf) : QbitBase({rf},cf) {}

  virtual ~QbitBase() {}

};


namespace qbit {

typedef std::shared_ptr<const QbitBase> Ptr;

inline const Tridiagonal sigmaop(Ptr qbit) {return mode::aop(qbit);}

const Tridiagonal sigmadagsigmaop();

const Tridiagonal sigmaxop(Ptr);
const Tridiagonal sigmayop(Ptr);
const Tridiagonal sigmazop();


inline double saturation(const StateVectorLow& psi) {return mode::photonNumber(psi);}

inline double saturation(const LazyDensityOperator& m) {return mode::photonNumber(m);}


inline StateVector state0() {return mode::fock(0,2);}
inline StateVector state1() {return mode::fock(1,2);}
StateVector init(dcomp psi1);
inline StateVector init(const Pars& p) {return init(p.qbitInit);}


Ptr make(const ParsDrivenDissipative&, QM_Picture);

Ptr make(const ParsDrivenDissipativePhaseNoise&, QM_Picture);



class Exact : public mode::Exact
{
public:
  Exact(dcomp zI) : mode::Exact(zI,0.,2) {}

};



template<bool IS_TIME_DEPENDENT>
class Hamiltonian : public mode::Hamiltonian<IS_TIME_DEPENDENT>
{
public:
  Hamiltonian(dcomp zSch, dcomp zI, dcomp eta)
    : mode::Hamiltonian<true >(zSch,zI,-eta,0,0,2) {}

  Hamiltonian(dcomp zSch,                  dcomp eta)
    : mode::Hamiltonian<false>(zSch,   -eta,0,0,2) {}

};



class Liouvillian : public mode::Liouvillian<false>
{
protected:
  Liouvillian(double gamma) : mode::Liouvillian<false>(gamma,0,keyTitle) {}

};


class LiouvillianPhaseNoise : public structure::ElementLiouvillianStrategies<1,2>
{
protected:
  LiouvillianPhaseNoise(double gamma_perpendicular, double gamma_parallel);
  
};

class LiouvillianDrivenPhaseNoise : public structure::ElementLiouvillianStrategies<1,3>
{
protected:
  LiouvillianDrivenPhaseNoise(double gamma_perpendicular, double gamma_pump, double gamma_parallel);
  
};


} // qbit



////////////////
//
// Highest level
//
////////////////


class Qbit
  : public qbit::Exact, public QbitBase
{
public:
  Qbit(const qbit::Pars&);
};

typedef Qbit QbitUIP;


class QbitSch
  : public qbit::Hamiltonian<false>, public QbitBase
{
public:
  QbitSch(const qbit::Pars&);
};


class DrivenQbit
  : public qbit::Hamiltonian<true>, public QbitBase
{
public:
  DrivenQbit(const qbit::ParsDriven&);
};

typedef DrivenQbit DrivenQbitUIP;


class DrivenQbitSch
  : public qbit::Hamiltonian<false>, public QbitBase
{
public:
  DrivenQbitSch(const qbit::ParsDriven&);
};


class DissipativeQbit
  : public qbit::Liouvillian, public qbit::Exact, public QbitBase
{
public:
  DissipativeQbit(double, double);

  template<typename BASE>
  DissipativeQbit(const qbit::ParsDissipative<BASE>& p) : DissipativeQbit(p.delta,p.gamma) {}
  
};


class DissipativeQbitSch
  : public qbit::Liouvillian, public qbit::Hamiltonian<false>, public QbitBase
{
public:
  DissipativeQbitSch(double, double);
  
  template<typename BASE>
  DissipativeQbitSch(const qbit::ParsDissipative<BASE>& p) : DissipativeQbitSch(p.delta,p.gamma) {}
};


class DissipativeQbitUIP
  : public qbit::Liouvillian, public qbit::Hamiltonian<true>, public QbitBase
{
public:
  DissipativeQbitUIP(double, double);

  template<typename BASE>
  DissipativeQbitUIP(const qbit::ParsDissipative<BASE>& p) : DissipativeQbitUIP(p.delta,p.gamma) {}
};


class DrivenDissipativeQbit 
  : public qbit::Liouvillian, public qbit::Hamiltonian<true>, public QbitBase
{
public:
  DrivenDissipativeQbit(const qbit::ParsDrivenDissipative& p);
};


class DrivenDissipativeQbitUIP
  : public qbit::Liouvillian, public qbit::Hamiltonian<true>, public QbitBase
{
public:
  DrivenDissipativeQbitUIP(const qbit::ParsDrivenDissipative&);
};


class DrivenDissipativeQbitSch
  : public qbit::Liouvillian, public qbit::Hamiltonian<false>, public QbitBase
{
public:
  DrivenDissipativeQbitSch(const qbit::ParsDrivenDissipative&);
};


class DissipativeQbitWithPhaseNoise
  : public qbit::Exact, public qbit::LiouvillianPhaseNoise, public QbitBase
{
public:
  DissipativeQbitWithPhaseNoise(double, double, double);

  template<typename BASE>
  DissipativeQbitWithPhaseNoise(const qbit::ParsDissipativePhaseNoise<BASE>& p) : DissipativeQbitWithPhaseNoise(p.delta,p.gamma,p.gamma_parallel) {}
  
};


class DissipativeQbitWithPhaseNoiseUIP
  : public qbit::Hamiltonian<true>, public qbit::LiouvillianPhaseNoise, public QbitBase
{
public:
  DissipativeQbitWithPhaseNoiseUIP(double, double, double);

  template<typename BASE>
  DissipativeQbitWithPhaseNoiseUIP(const qbit::ParsDissipativePhaseNoise<BASE>& p) : DissipativeQbitWithPhaseNoiseUIP(p.delta,p.gamma,p.gamma_parallel) {}
  
};

class DrivenDissipativeQbitWithPhaseNoise
  : public qbit::Hamiltonian<true>, public qbit::LiouvillianPhaseNoise, public QbitBase
{
public:
  DrivenDissipativeQbitWithPhaseNoise(const qbit::ParsDrivenDissipativePhaseNoise&);
};

class DrivenDissipativeQbitWithPhaseNoiseUIP
  : public qbit::Hamiltonian<true>, public qbit::LiouvillianPhaseNoise, public QbitBase
{
public:
  DrivenDissipativeQbitWithPhaseNoiseUIP(const qbit::ParsDrivenDissipativePhaseNoise&);
};



class DissipativeQbitWithIncoherentPumpAndPhaseNoise
  : public qbit::Exact, public qbit::LiouvillianDrivenPhaseNoise, public QbitBase
{
public:
  DissipativeQbitWithIncoherentPumpAndPhaseNoise(double, double, double, double);

  template<typename PARS>
  DissipativeQbitWithIncoherentPumpAndPhaseNoise(PARS p) : DissipativeQbitWithIncoherentPumpAndPhaseNoise(p.delta,p.gamma,p.gamma_pump,p.gamma_parallel) {}
  
};


class DissipativeQbitWithIncoherentPumpAndPhaseNoiseUIP
  : public qbit::Hamiltonian<true>, public qbit::LiouvillianDrivenPhaseNoise, public QbitBase
{
public:
  DissipativeQbitWithIncoherentPumpAndPhaseNoiseUIP(double, double, double, double);

  template<typename PARS>
  DissipativeQbitWithIncoherentPumpAndPhaseNoiseUIP(PARS p) : DissipativeQbitWithIncoherentPumpAndPhaseNoiseUIP(p.delta,p.gamma,p.gamma_pump,p.gamma_parallel) {}
  
};


#endif // CPPQEDELEMENTS_FREES_QBIT__H_INCLUDED
