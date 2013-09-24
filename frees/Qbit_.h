// -*- C++ -*-
#ifndef ELEMENTS_FREES_QBIT__H_INCLUDED
#define ELEMENTS_FREES_QBIT__H_INCLUDED

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
  void           process_v(        Averages&                 ) const {}

};


} // qbit




class QbitBase 
  : public ModeBase, public qbit::Averaged
{
public:
  QbitBase(const RealFreqs& realFreqs=RealFreqs(), const ComplexFreqs& complexFreqs=ComplexFreqs());

  virtual ~QbitBase() {}

};

namespace qbit {

using namespace structure::free;

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


#endif // ELEMENTS_FREES_QBIT__H_INCLUDED
