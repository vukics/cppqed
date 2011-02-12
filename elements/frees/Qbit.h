// -*- C++ -*-
#ifndef _ELEMENT_QBIT_H
#define _ELEMENT_QBIT_H

#include "QbitFwd.h"

#include "ParsQbit.h"

#include "QM_PictureFwd.h"
#include "StateVectorFwd.h"

#include "ElementLiouvillean.h"
#include "ElementAveraged.h"
#include "Free.h"
#include "FreeExact.h"
#include "TridiagonalHamiltonian.h"

#include<boost/shared_ptr.hpp>


namespace qbit {

using namespace structure::free;

typedef boost::shared_ptr<const QbitBase> SmartPtr;

const Tridiagonal sigmaop();
const Tridiagonal sigmadagop();
const Tridiagonal sigmadagsigmaop();
const Tridiagonal sigmaxop();
const Tridiagonal sigmayop();
const Tridiagonal sigmazop();

const Frequencies freqs(const dcomp& zI);

const Frequencies freqs(const QbitBase*);


double saturation(const StateVectorLow&); 

double saturation(const LazyDensityOperator&);


const StateVector state0();
const StateVector state1();
const StateVector init(const dcomp& psi1);
const StateVector init(const Pars&);


SmartPtr maker(const ParsPumpedLossy&, QM_Picture);




class Exact : public structure::FreeExact
{
public:
  Exact(const dcomp& zI);

  const dcomp& get_zI() const {return zI_;}
  
private:
  void updateU(double) const;

  bool isUnitary() const {return !bool(real(zI_));}

  const dcomp zI_;

};



namespace details {

struct Empty {};

}


template<bool IS_TD>
class Hamiltonian 
  : public structure::TridiagonalHamiltonian<1,IS_TD>,
    public mpl::if_c<IS_TD,Exact,details::Empty>::type
{
public:
  typedef structure::TridiagonalHamiltonian<1,IS_TD> Base;

  Hamiltonian(const dcomp& zSch, const dcomp& zI, const dcomp& eta, mpl::bool_<IS_TD> =mpl:: true_()); // works for IS_TD=true
  Hamiltonian(const dcomp& zSch,                  const dcomp& eta, mpl::bool_<IS_TD> =mpl::false_()); // works for IS_TD=false

protected:
  const dcomp& get_zSch() const {return zSch_;}
  const dcomp& get_zI  () const {return Exact::get_zI();}

  const dcomp& get_eta() const {return eta_;}

private:
  const dcomp zSch_, eta_;

};



class Liouvillean : public structure::ElementLiouvillean<1,1>
{
protected:
  Liouvillean(double gamma) : gamma_(gamma) {}

private:
  void   doActWithJ (StateVectorLow&           ) const;
  double probability(const LazyDensityOperator&) const;

  const double gamma_;

};
 


class Averaged
  : public structure::ElementAveraged<1>
{
public:
  typedef structure::ElementAveraged<1> Base;

  Averaged();

private:
  const Averages average(const LazyDensityOperator&) const;
  void           process(Averages&                 ) const {}

};


} // qbit



////////////////
//
// Highest level
//
////////////////

class QbitBase 
  : public structure::Free, public qbit::Averaged
{
public:
  QbitBase(const RealFreqs& realFreqs=RealFreqs(), const ComplexFreqs& complexFreqs=ComplexFreqs());

  virtual ~QbitBase() {}

};


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


#endif // _ELEMENT_QBIT_H
