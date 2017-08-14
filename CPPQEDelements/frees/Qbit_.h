// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines the qbit-bundle (tackling the dynamics of a single qbit)}
#ifndef CPPQEDELEMENTS_FREES_QBIT__H_INCLUDED
#define CPPQEDELEMENTS_FREES_QBIT__H_INCLUDED

#include "Qbit_Fwd.h"

#include "ParsQbit.h"

#include "Mode_.h"
#include "StateVector.h"

// Qbit traced back to Mode with dim=2

namespace qbit {

using mode::NoTime;

const std::string keyTitle="Qbit";


class Averaged
  : public structure::ElementAveraged<1>
{
public:
  typedef structure::ElementAveraged<1> Base;

  Averaged();

private:
  const Averages average_v(NoTime, const LazyDensityOperator&) const;

};


} // qbit




class QbitBase 
  : public ModeBase, public qbit::Averaged
{
public:
  explicit QbitBase(const RealFreqs& =emptyRF, const ComplexFreqs& =emptyCF);

  explicit QbitBase(const ComplexFreqs& complexFreqs) : QbitBase(emptyRF,complexFreqs) {}
  explicit QbitBase(RealFreqsInitializer rf, ComplexFreqsInitializer cf={}) : QbitBase(RealFreqs(rf),ComplexFreqs(cf)) {}
  explicit QbitBase(ComplexFreqsInitializer cf) : QbitBase({},cf) {}
  explicit QbitBase(RF rf, CF cf=CF()) : QbitBase(RealFreqsInitializer{rf}, cf==CF() ? ComplexFreqsInitializer{} : ComplexFreqsInitializer{cf}) {}
  explicit QbitBase(CF cf) : QbitBase(ComplexFreqsInitializer{cf}) {}
  explicit QbitBase(RealFreqsInitializer rf, CF cf) : QbitBase(rf,{cf}) {}
  explicit QbitBase(RF rf, ComplexFreqsInitializer cf) : QbitBase({rf},cf) {}

  virtual ~QbitBase() {}

};


namespace qbit {

using namespace structure::freesystem;

typedef boost::shared_ptr<const QbitBase> Ptr;

inline const Tridiagonal sigmaop(Ptr qbit) {return mode::aop(qbit);}

const Tridiagonal sigmadagsigmaop();

const Tridiagonal sigmaxop(Ptr);
const Tridiagonal sigmayop(Ptr);
const Tridiagonal sigmazop();


inline double saturation(const StateVectorLow& psi) {return mode::photonNumber(psi);}

inline double saturation(const LazyDensityOperator& m) {return mode::photonNumber(m);}


inline const StateVector state0() {return mode::fock(0,2);}
inline const StateVector state1() {return mode::fock(1,2);}
const StateVector init(const dcomp& psi1);
inline const StateVector init(const Pars& p) {return init(p.qbitInit);}


Ptr make(const ParsPumpedLossy&, QM_Picture);




class Exact : public mode::Exact
{
public:
  Exact(const dcomp& zI) : mode::Exact(zI,2) {}

};



template<bool IS_TIME_DEPENDENT>
class Hamiltonian : public mode::Hamiltonian<IS_TIME_DEPENDENT>
{
public:
  Hamiltonian(const dcomp& zSch, const dcomp& zI, const dcomp& eta)
    : mode::Hamiltonian<true >(zSch,zI,-eta,2) {}

  Hamiltonian(const dcomp& zSch,                  const dcomp& eta)
    : mode::Hamiltonian<false>(zSch,   -eta,2) {}

};



class Liouvillean : public mode::Liouvillean<false>
{
protected:
  Liouvillean(double gamma) : mode::Liouvillean<false>(gamma,0,keyTitle) {}

};


class LiouvilleanPhaseNoise : public structure::ElementLiouvilleanStrategies<1,2>
{
protected:
  LiouvilleanPhaseNoise(double gamma_perpendicular, double gamma_parallel);
  
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


class PumpedQbit
  : public qbit::Hamiltonian<true>, public QbitBase
{
public:
  PumpedQbit(const qbit::ParsPumped&);
};

typedef PumpedQbit PumpedQbitUIP;


class PumpedQbitSch
  : public qbit::Hamiltonian<false>, public QbitBase
{
public:
  PumpedQbitSch(const qbit::ParsPumped&);
};


class LossyQbit
  : public qbit::Liouvillean, public qbit::Exact, public QbitBase
{
public:
  LossyQbit(const qbit::ParsLossy&);
};


class LossyQbitSch
  : public qbit::Liouvillean, public qbit::Hamiltonian<false>, public QbitBase
{
public:
  LossyQbitSch(const qbit::ParsLossy&);
};


class LossyQbitUIP
  : public qbit::Liouvillean, public qbit::Hamiltonian<true>, public QbitBase
{
public:
  LossyQbitUIP(const qbit::ParsLossy&);
};


class PumpedLossyQbit 
  : public qbit::Liouvillean, public qbit::Hamiltonian<true>, public QbitBase
{
public:
  PumpedLossyQbit(const qbit::ParsPumpedLossy&);
};


class PumpedLossyQbitUIP
  : public qbit::Liouvillean, public qbit::Hamiltonian<true>, public QbitBase
{
public:
  PumpedLossyQbitUIP(const qbit::ParsPumpedLossy&);
};


class PumpedLossyQbitSch
  : public qbit::Liouvillean, public qbit::Hamiltonian<false>, public QbitBase
{
public:
  PumpedLossyQbitSch(const qbit::ParsPumpedLossy&);
};


class LossyQbitWithPhaseNoise
  : public qbit::Exact, public qbit::LiouvilleanPhaseNoise, public QbitBase
{
public:
  LossyQbitWithPhaseNoise(const qbit::ParsLossy&, double gamma_parallel);
  
};


class LossyQbitWithPhaseNoiseUIP
  : public qbit::Hamiltonian<true>, public qbit::LiouvilleanPhaseNoise, public QbitBase
{
public:
  LossyQbitWithPhaseNoiseUIP(const qbit::ParsLossy&, double gamma_parallel);
  
};

class PumpedLossyQbitWithPhaseNoise
  : public qbit::Hamiltonian<true>, public qbit::LiouvilleanPhaseNoise, public QbitBase
{
public:
  PumpedLossyQbitWithPhaseNoise(const qbit::ParsPumpedLossy&, double gamma_parallel);
};

class PumpedLossyQbitWithPhaseNoiseUIP
  : public qbit::Hamiltonian<true>, public qbit::LiouvilleanPhaseNoise, public QbitBase
{
public:
  PumpedLossyQbitWithPhaseNoiseUIP(const qbit::ParsPumpedLossy&, double gamma_parallel);
};

#endif // CPPQEDELEMENTS_FREES_QBIT__H_INCLUDED
