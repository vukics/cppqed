// -*- C++ -*-
#ifndef _MODE__H
#define _MODE__H

#include "Mode_Fwd.h"

#include "ParsMode.h" // for the UI

#include "QM_PictureFwd.h"
#include "StateVectorFwd.h"

#include "ElementLiouvillean.h"
#include "ElementAveraged.h"
#include "Free.h"
#include "FreeExact.h"
#include "TridiagonalHamiltonian.h"

#include<boost/shared_ptr.hpp>

namespace mode {

const std::string keyTitle="Mode";


using namespace structure::free;

typedef boost::shared_ptr<const ModeBase> SmartPtr;

const Tridiagonal aop(const ModeBase*);
const Tridiagonal nop(const ModeBase*);

inline const Tridiagonal xop(const ModeBase* mode) {return tridiagPlusHC(aop(mode))/sqrt(2.);}
// inline const Tridiagonal yop(const ModeBase*) {return ...}

struct PrepError : public cpputils::Exception {};

const StateVector coherent(const dcomp&, size_t);
const StateVector fock(size_t n, size_t dim, double phase=0) throw(PrepError);
const StateVector init(const Pars&);


template<typename A>
const SmartPtr make(const Pars           &, QM_Picture, const A&);

template<typename A>
const SmartPtr make(const ParsLossy      &, QM_Picture, const A&);

template<typename A>
const SmartPtr make(const ParsPumped     &, QM_Picture, const A&);

template<typename A>
const SmartPtr make(const ParsPumpedLossy&, QM_Picture, const A&);


double photonNumber(const StateVectorLow&); 
// This can be implemented in an extremely elegant way using tensor
// notation from blitz++, but it is not enough, the one below is
// needed.
double photonNumber(const LazyDensityOperator&);


inline std::ostream& isFiniteTempStream(std::ostream& os, double    , boost::mpl::false_)
{return os                                     <<std::endl;}
inline std::ostream& isFiniteTempStream(std::ostream& os, double nTh, boost::mpl:: true_)
{return os<<" Finite temperature.\n# nTh="<<nTh<<std::endl;}


inline double finiteTemperatureHamiltonianDecay(const ParsLossy& p, boost::mpl::false_) {return p.kappa             ;}
inline double finiteTemperatureHamiltonianDecay(const ParsLossy& p, boost::mpl:: true_) {return p.kappa*(2.*p.nTh+1);}


const Tridiagonal::Diagonal mainDiagonal(const dcomp& z, size_t dim);

const Tridiagonal pumping(const dcomp& eta, size_t dim);



class Exact : public structure::FreeExact
{
public:
  Exact(const dcomp& zI, size_t);

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

  Hamiltonian(const dcomp& zSch, const dcomp& zI, const dcomp& eta, size_t, mpl::bool_<IS_TD> =mpl:: true_()); // works for IS_TD=true
  Hamiltonian(const dcomp& zSch,                  const dcomp& eta, size_t, mpl::bool_<IS_TD> =mpl::false_()); // works for IS_TD=false
  // The trailing dummy argument is there to cause a _compile_time_
  // (and not merely linking time) error in case of misuse

protected:
  const dcomp& get_zSch() const {return zSch_;}
  const dcomp& get_zI  () const {return Exact::get_zI();}

  const dcomp& get_eta() const {return eta_;}

private:
  const dcomp zSch_,
  // z_=kappa_-I*delta_ --- zSch_ is the part appearing in the
  // tridiagonal hamiltonian, zI_ is appearing in the frequencies.
    eta_;
  const size_t dim_;

};



template<bool IS_FINITE_TEMP, bool IS_ALTERNATIVE=false> class Liouvillean;


template<bool IS_ALTERNATIVE> 
class Liouvillean<false,IS_ALTERNATIVE> 
  : protected boost::mpl::false_, // Tagging for the isFiniteTemp... functions
    public structure::ElementLiouvillean<1,1>
{
protected:
  Liouvillean(double kappa, double=0, const std::string& kT=keyTitle) : structure::ElementLiouvillean<1,1>(kT,"excitation loss"), kappa_(kappa) {}
  // the trailing dummy argument is there only to have the same form for the ctor as in the IS_FINITE_TEMP=true case

private:
  void   doActWithJ (StateVectorLow&           ) const;
  double probability(const LazyDensityOperator&) const;

  const double kappa_;

};
 

template<> 
class Liouvillean<true >
  : protected boost::mpl::true_,
    public structure::ElementLiouvillean<1,2>
{
protected:
  typedef structure::ElementLiouvillean<1,2> Base;

  Liouvillean(double kappa, double nTh, const std::string& kT=keyTitle);

private:
  const double kappa_, nTh_;

};



// A basic, extensible Averaging class
class Averaged : public structure::ClonableElementAveraged<1>
{
public:
  typedef structure::ClonableElementAveraged<1> Base;

  Averaged(const KeyLabels& follow=KeyLabels(), const KeyLabels& precede=KeyLabels());

  virtual const Averages average(const LazyDensityOperator&) const;
  virtual void           process(Averages&)                  const;

private:
  const ClonedPtr do_clone() const {return new Averaged(*this);}

};


inline const SmartPtr make(const Pars           & p, QM_Picture qmp) {return make(p,qmp,Averaged());}
inline const SmartPtr make(const ParsLossy      & p, QM_Picture qmp) {return make(p,qmp,Averaged());}
inline const SmartPtr make(const ParsPumped     & p, QM_Picture qmp) {return make(p,qmp,Averaged());}
inline const SmartPtr make(const ParsPumpedLossy& p, QM_Picture qmp) {return make(p,qmp,Averaged());}



class AveragedQuadratures : public Averaged
{
public:
  AveragedQuadratures(const KeyLabels& follow=KeyLabels(), const KeyLabels& precede=KeyLabels());

  virtual const Averages average(const LazyDensityOperator&) const;
  virtual void           process(Averages&)                  const;

private:
  const ClonedPtr do_clone() const {return new AveragedQuadratures(*this);}

};


template<typename Base>
class AveragedMonitorCutoff : public Base
{
public:
  typedef typename Base::Averages Averages;
  typedef typename Base::KeyLabels KeyLabels;

  AveragedMonitorCutoff();

  const Averages average(const LazyDensityOperator&) const;
  void           process(Averages&)                  const;

};


} // mode



////////////////
//
// Highest level
//
////////////////

// Notice how everything here is achieved with trivial class composition


class ModeBase 
  : public structure::Free
{
public:
  ModeBase(size_t dim,
	   const RealFreqs& realFreqs=RealFreqs(), const ComplexFreqs& complexFreqs=ComplexFreqs(),
	   const std::string& keyTitle=mode::keyTitle);

  virtual ~ModeBase() {}

};


/////
// 1
/////

template<typename A>
class Mode 
  : public mode::Exact, public ModeBase, public A
{
public:
  Mode(const mode::Pars&, const A& =A());
};


template<typename A>
struct ModeUIP : Mode<A>
// in this case the uip and ip coincide, 
{
  ModeUIP(const mode::Pars& p, const A& a=A()) : Mode<A>(p,a) {}
};


/////
// 2
/////


template<typename A>
class ModeSch
  : public mode::Hamiltonian<false>, public ModeBase, public A
{
public:
  ModeSch(const mode::Pars&, const A& =A());
};


/////
// 3
/////

template<typename A>
class PumpedMode 
  : public mode::Hamiltonian<true>, public ModeBase, public A
{
public:
  PumpedMode(const mode::ParsPumped&, const A& =A());
};


template<typename A>
struct PumpedModeUIP : PumpedMode<A>
{
  PumpedModeUIP(const mode::ParsPumped& p, const A& a=A()) : PumpedMode<A>(p,a) {}
};


/////
// 4
/////

template<typename A>
class PumpedModeSch
  : public mode::Hamiltonian<false>, public ModeBase, public A
{
public:
  PumpedModeSch(const mode::ParsPumped&, const A& =A());
};


/////
// 5
/////

template<bool IS_FINITE_TEMP, typename A>
class LossyMode 
  : public mode::Liouvillean<IS_FINITE_TEMP>, public mode::Exact, public ModeBase, public A
{
public:
  LossyMode(const mode::ParsLossy&, const A& =A());

};

/////
// 6
/////

template<bool IS_FINITE_TEMP, typename A>
class LossyModeUIP 
  : public mode::Liouvillean<IS_FINITE_TEMP>, public mode::Hamiltonian<true>, public ModeBase, public A
{
public:
  LossyModeUIP(const mode::ParsLossy&, const A& =A());

};

/////
// 7
/////

template<bool IS_FINITE_TEMP, typename A>
class LossyModeSch 
  : public mode::Liouvillean<IS_FINITE_TEMP>, public mode::Hamiltonian<false>, public ModeBase, public A
{
public:
  LossyModeSch(const mode::ParsLossy&, const A& =A());

};


/////
// 8
/////

template<bool IS_FINITE_TEMP, typename A>
class PumpedLossyMode 
  : public mode::Liouvillean<IS_FINITE_TEMP>, public mode::Hamiltonian<true>, public ModeBase, public A
{
public:
  PumpedLossyMode(const mode::ParsPumpedLossy&, const A& =A());
};

/////
// 9
/////

template<bool IS_FINITE_TEMP, typename A>
class PumpedLossyModeUIP 
  : public mode::Liouvillean<IS_FINITE_TEMP>, public mode::Hamiltonian<true>, public ModeBase, public A
{
public:
  PumpedLossyModeUIP(const mode::ParsPumpedLossy&, const A& =A());
};

/////
// 10
/////

template<bool IS_FINITE_TEMP, typename A>
class PumpedLossyModeSch
  : public mode::Liouvillean<IS_FINITE_TEMP>, public mode::Hamiltonian<false>, public ModeBase, public A
{
public:
  PumpedLossyModeSch(const mode::ParsPumpedLossy&, const A& =A());
};


//////////////////////////////////////////////////////////////////////
// One more to help testing the alternative jumping in MCWF_Trajectory
//////////////////////////////////////////////////////////////////////

template<typename A=mode::Averaged>
class PumpedLossyModeAlternative 
  : public mode::Liouvillean<false,true>, public mode::Hamiltonian<true>, public ModeBase, public A
{
public:
  PumpedLossyModeAlternative(const mode::ParsPumpedLossy&, const A& =A());
};


//////////////////////////////////////////////////////////////////////
// One more to test time-dependent jump and averages 
//////////////////////////////////////////////////////////////////////

class PumpedLossyModeIP_NoExact
  : public ModeBase,
    public structure::TridiagonalHamiltonian<1,true>, 
    public structure::ElementLiouvillean<1,1,true>,
    public structure::ElementAveraged<1,true>
{
public:
  PumpedLossyModeIP_NoExact(const mode::ParsPumpedLossy&);

private:
  typedef structure::TridiagonalHamiltonian<1,true>::StateVectorLow StateVectorLow;
  typedef structure::ElementLiouvillean<1,1,true>::LazyDensityOperator LazyDensityOperator;

  void   doActWithJ (double t, StateVectorLow&           ) const;
  double probability(double t, const LazyDensityOperator&) const;

  const Averages average(double t, const LazyDensityOperator&) const;

  void           process(Averages&)                  const {}

  const dcomp z_;

};


#endif // _MODE__H
