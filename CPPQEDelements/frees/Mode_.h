// -*- C++ -*-
#ifndef FREES_MODE__H_INCLUDED
#define FREES_MODE__H_INCLUDED

#include "Mode_Fwd.h"

#include "ParsMode.h" // for the user interface

#include "QM_PictureFwd.h"
#include "StateVectorFwd.h"

#include "ElementLiouvillean.h"
#include "ElementAveraged.h"
#include "Free.h"
#include "FreeExact.h"
#include "TridiagonalHamiltonian.h"

#include <boost/shared_ptr.hpp>

namespace mode {

const std::string keyTitle="Mode";


using namespace structure::free; using structure::NoTime;

typedef boost::shared_ptr<const ModeBase> Ptr;

const Tridiagonal aop(Ptr);
const Tridiagonal nop(Ptr);

inline const Tridiagonal xop(Ptr mode) {return tridiagPlusHC(aop(mode))/sqrt(2.);}
// inline const Tridiagonal yop(mode::Ptr) {return ...}

struct PrepError : public cpputils::Exception {};

const StateVector coherent(const dcomp&, size_t);
const StateVector fock(size_t n, size_t dim, double phase=0) throw(PrepError);
const StateVector init(const Pars&);


template<typename AveragingType, typename... AveragingConstructorParameters>
const Ptr make(const Pars           &, QM_Picture, const AveragingConstructorParameters&... );

template<typename AveragingType, typename... AveragingConstructorParameters>
const Ptr make(const ParsLossy      &, QM_Picture, const AveragingConstructorParameters&... );

template<typename AveragingType, typename... AveragingConstructorParameters>
const Ptr make(const ParsPumped     &, QM_Picture, const AveragingConstructorParameters&... );

template<typename AveragingType, typename... AveragingConstructorParameters>
const Ptr make(const ParsPumpedLossy&, QM_Picture, const AveragingConstructorParameters&... );


const Ptr make(const Pars           &, QM_Picture);
const Ptr make(const ParsLossy      &, QM_Picture);
const Ptr make(const ParsPumped     &, QM_Picture);
const Ptr make(const ParsPumpedLossy&, QM_Picture);


double photonNumber(const StateVectorLow&); 
// This can be implemented in an extremely elegant way using tensor notation from blitz++, but it is not enough, the one below is needed.
double photonNumber(const LazyDensityOperator&);


inline std::ostream& isFiniteTempStream(std::ostream& os, double    , boost::mpl::false_)
{return os                                     <<std::endl;}
inline std::ostream& isFiniteTempStream(std::ostream& os, double nTh, boost::mpl:: true_)
{return os<<" Finite temperature.\n# nTh="<<nTh<<std::endl;}


inline double finiteTemperatureHamiltonianDecay(const ParsLossy& p, boost::mpl::false_) {return p.kappa             ;}
inline double finiteTemperatureHamiltonianDecay(const ParsLossy& p, boost::mpl:: true_) {return p.kappa*(2.*p.nTh+1);}


const Tridiagonal::Diagonal mainDiagonal(const dcomp& z, size_t dim);

const Tridiagonal pumping(const dcomp& eta, size_t dim);



class Exact : public structure::FreeExact<false>
{
public:
  Exact(const dcomp& zI, size_t);

  const dcomp& get_zI() const {return zI_;}
  
private:
  void updateU(Time) const; ///< `Time` is structure::OneTime in this case

  bool isUnitary_v() const {return !bool(real(zI_));}

  const dcomp zI_;

};



template<bool IS_TIME_DEPENDENT>
class Hamiltonian 
  : public structure::TridiagonalHamiltonian<1,IS_TIME_DEPENDENT>,
    public mpl::if_c<IS_TIME_DEPENDENT,Exact,mpl::empty_base>::type
{
public:
  typedef structure::TridiagonalHamiltonian<1,IS_TIME_DEPENDENT> Base;

  Hamiltonian(const dcomp& zSch, const dcomp& zI, const dcomp& eta, size_t, mpl::bool_<IS_TIME_DEPENDENT> =mpl:: true_()); // works for IS_TIME_DEPENDENT=true
  Hamiltonian(const dcomp& zSch,                  const dcomp& eta, size_t, mpl::bool_<IS_TIME_DEPENDENT> =mpl::false_()); // works for IS_TIME_DEPENDENT=false
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
  void   doActWithJ (NoTime, StateVectorLow&           ) const;
  double rate       (NoTime, const LazyDensityOperator&) const;

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
  void doActWithJ(NoTime, StateVectorLow&, JumpNo<0>) const;
  void doActWithJ(NoTime, StateVectorLow&, JumpNo<1>) const;
  
  double rate(NoTime, const LazyDensityOperator&, JumpNo<0>) const;
  double rate(NoTime, const LazyDensityOperator&, JumpNo<1>) const;
  
  const double kappa_, nTh_;

};



// A basic, extensible Averaging class
class Averaged : public structure::ClonableElementAveraged<1>
{
public:
  typedef structure::ClonableElementAveraged<1> Base;

  Averaged(const KeyLabels& follow=KeyLabels(), const KeyLabels& precede=KeyLabels());

protected:
  const Averages average_v(NoTime, const LazyDensityOperator&) const;
  void           process_v(        Averages&)                  const;

private:
  const ClonedPtr do_clone() const {return new Averaged(*this);}

};


class AveragedQuadratures : public Averaged
{
public:
  AveragedQuadratures(const KeyLabels& follow=KeyLabels(), const KeyLabels& precede=KeyLabels());

protected:
  const Averages average_v(NoTime, const LazyDensityOperator&) const;
  void           process_v(        Averages&)                  const;

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

private:
  const Averages average_v(NoTime, const LazyDensityOperator&) const;
  void           process_v(        Averages&)                  const;

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
           const RealFreqs& realFreqs=emptyRF, const ComplexFreqs& complexFreqs=emptyCF,
           const std::string& keyTitle=mode::keyTitle);

  ModeBase(size_t dim, const ComplexFreqs& complexFreqs, const std::string& keyTitle=mode::keyTitle) : ModeBase(dim,emptyRF,complexFreqs,keyTitle) {}
  ModeBase(size_t dim, RealFreqsInitializer rf, ComplexFreqsInitializer cf={}, const std::string& keyTitle=mode::keyTitle)
    : ModeBase(dim,RealFreqs(rf),ComplexFreqs(cf),keyTitle) {}
  ModeBase(size_t dim, ComplexFreqsInitializer cf, const std::string& keyTitle=mode::keyTitle) : ModeBase(dim,{},cf,keyTitle) {}
  ModeBase(size_t dim, RF rf, CF cf=CF(), const std::string& keyTitle=mode::keyTitle)
    : ModeBase(dim,RealFreqsInitializer{rf}, cf==CF() ? ComplexFreqsInitializer{} : ComplexFreqsInitializer{cf},keyTitle) {}
  ModeBase(size_t dim, CF cf, const std::string& keyTitle=mode::keyTitle) : ModeBase(dim,ComplexFreqsInitializer{cf},keyTitle) {}
  ModeBase(size_t dim, RealFreqsInitializer rf, CF cf, const std::string& keyTitle=mode::keyTitle) : ModeBase(dim,rf,{cf},keyTitle) {}
  ModeBase(size_t dim, RF rf, ComplexFreqsInitializer cf, const std::string& keyTitle=mode::keyTitle) : ModeBase(dim,{rf},cf,keyTitle) {}

  virtual ~ModeBase() {}

};


/////
// 1
/////

template<typename AveragingType>
class Mode 
  : public mode::Exact, public ModeBase, public AveragingType
{
public:
  template<typename... AveragingConstructorParameters>
  Mode(const mode::Pars&, const AveragingConstructorParameters&... );
};


template<typename AveragingType>
struct ModeUIP : Mode<AveragingType>
// in this case the uip and ip coincide, 
{
  template<typename... AveragingConstructorParameters>
  ModeUIP(const mode::Pars& p, const AveragingConstructorParameters&... a) : Mode<AveragingType>(p,a...) {}
};


/////
// 2
/////


template<typename AveragingType>
class ModeSch
  : public mode::Hamiltonian<false>, public ModeBase, public AveragingType
{
public:
  template<typename... AveragingConstructorParameters>
  ModeSch(const mode::Pars&, const AveragingConstructorParameters&... );
};


/////
// 3
/////

template<typename AveragingType>
class PumpedMode 
  : public mode::Hamiltonian<true>, public ModeBase, public AveragingType
{
public:
  template<typename... AveragingConstructorParameters>
  PumpedMode(const mode::ParsPumped&, const AveragingConstructorParameters&... );
};


template<typename AveragingType>
struct PumpedModeUIP : PumpedMode<AveragingType>
{
  template<typename... AveragingConstructorParameters>
  PumpedModeUIP(const mode::ParsPumped& p, const AveragingConstructorParameters&... a) : PumpedMode<AveragingType>(p,a...) {}
};


/////
// 4
/////

template<typename AveragingType>
class PumpedModeSch
  : public mode::Hamiltonian<false>, public ModeBase, public AveragingType
{
public:
  template<typename... AveragingConstructorParameters>
  PumpedModeSch(const mode::ParsPumped&, const AveragingConstructorParameters&... );
};


/////
// 5
/////

template<bool IS_FINITE_TEMP, typename AveragingType>
class LossyMode 
  : public mode::Liouvillean<IS_FINITE_TEMP>, public mode::Exact, public ModeBase, public AveragingType
{
public:
  template<typename... AveragingConstructorParameters>
  LossyMode(const mode::ParsLossy&, const AveragingConstructorParameters&... );

};

/////
// 6
/////

template<bool IS_FINITE_TEMP, typename AveragingType>
class LossyModeUIP 
  : public mode::Liouvillean<IS_FINITE_TEMP>, public mode::Hamiltonian<true>, public ModeBase, public AveragingType
{
public:
  template<typename... AveragingConstructorParameters>
  LossyModeUIP(const mode::ParsLossy&, const AveragingConstructorParameters&... );

};

/////
// 7
/////

template<bool IS_FINITE_TEMP, typename AveragingType>
class LossyModeSch 
  : public mode::Liouvillean<IS_FINITE_TEMP>, public mode::Hamiltonian<false>, public ModeBase, public AveragingType
{
public:
  template<typename... AveragingConstructorParameters>
  LossyModeSch(const mode::ParsLossy&, const AveragingConstructorParameters&... );

};


/////
// 8
/////

template<bool IS_FINITE_TEMP, typename AveragingType>
class PumpedLossyMode 
  : public mode::Liouvillean<IS_FINITE_TEMP>, public mode::Hamiltonian<true>, public ModeBase, public AveragingType
{
public:
  template<typename... AveragingConstructorParameters>
  PumpedLossyMode(const mode::ParsPumpedLossy&, const AveragingConstructorParameters&... );
};

/////
// 9
/////

template<bool IS_FINITE_TEMP, typename AveragingType>
class PumpedLossyModeUIP 
  : public mode::Liouvillean<IS_FINITE_TEMP>, public mode::Hamiltonian<true>, public ModeBase, public AveragingType
{
public:
  template<typename... AveragingConstructorParameters>
  PumpedLossyModeUIP(const mode::ParsPumpedLossy&, const AveragingConstructorParameters&... );
};

/////
// 10
/////

template<bool IS_FINITE_TEMP, typename AveragingType>
class PumpedLossyModeSch
  : public mode::Liouvillean<IS_FINITE_TEMP>, public mode::Hamiltonian<false>, public ModeBase, public AveragingType
{
public:
  template<typename... AveragingConstructorParameters>
  PumpedLossyModeSch(const mode::ParsPumpedLossy&, const AveragingConstructorParameters&... );
};


//////////////////////////////////////////////////////////////////////
// One more to help testing the alternative jumping in MCWF_Trajectory
//////////////////////////////////////////////////////////////////////

template<typename AveragingType=mode::Averaged>
class PumpedLossyModeAlternative 
  : public mode::Liouvillean<false,true>, public mode::Hamiltonian<true>, public ModeBase, public AveragingType
{
public:
  template<typename... AveragingConstructorParameters>
  PumpedLossyModeAlternative(const mode::ParsPumpedLossy&, const AveragingConstructorParameters&... );
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
  typedef structure::OneTime OneTime;
  typedef structure::TridiagonalHamiltonian<1,true>::StateVectorLow StateVectorLow;
  typedef structure::ElementLiouvillean<1,1,true>::LazyDensityOperator LazyDensityOperator;

  void   doActWithJ (OneTime, StateVectorLow&           ) const;
  double rate       (OneTime, const LazyDensityOperator&) const;

  const Averages average_v(OneTime, const LazyDensityOperator&) const;

  const dcomp z_;

};

#endif // FREES_MODE__H_INCLUDED
