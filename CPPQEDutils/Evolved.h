// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_UTILS_EVOLVED_H_INCLUDED
#define CPPQEDCORE_UTILS_EVOLVED_H_INCLUDED

#include "MathExtensions.h"

#ifdef BZ_HAVE_BOOST_SERIALIZATION
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>
#endif // BZ_HAVE_BOOST_SERIALIZATION

#include <functional>
#include <iostream>
#include <memory>

/// Comprises utilities related to ODE adaptive evolution
namespace evolved {


/// Enumeration for different stepping-function types, for i/o operations
/**
 * At the moment, Runge-Kutta Cash-Karp (4,5) and Runge-Kutta Prince-Dormand (8,9) are defined, cf. http://www.gnu.org/software/gsl/manual/html_node/Stepping-Functions.html
 */
enum SteppingFunction {SF_RKCK, SF_RK8PD};

/// \name I/O operations for SteppingFunction
//@{
std::ostream& operator<<(std::ostream&, SteppingFunction);
std::istream& operator>>(std::istream&, SteppingFunction&);
//@}


/// Bookkeeps the timestep-data of Evolved
/**
 * Very rarely used in itself, for the most part it can be considered as a template-parameter independent base of Evolved.
 *
 * The timestep-regulation policy of adaptive-stepsize drivers is that for the step in hand they perform a timestep of *at most* `dtTry`. Thereupon, they report the performed timestep
 * (this is `dtDid`), and suggest a stepsize *for the next step to try*, which can be larger than *dtDid* (this will be the new value of `dtTry` after the step). Hence, the step can
 * diminish in the actual step, but can grow only in the next step.
 *
 * Timestep is regulated by two parameters, `epsRel` limiting the (element-wise) relative precision, and `epsAbs` the absolute precision. In general, it can be said that anything below
 * `epsAbs` in absolute value should be considered numerical trash.
 *
 * ODE steppers must be supplied with a sensible value for the initial timestep to try in the first step.
 *
 */
class TimeStepBookkeeper
{
public:
  /// \name Actual time getter/setter
  //@{
  double getTime(        ) const {return t_;}
  void   setTime(double t)       {t_=t;}
  //@}

  double getDtDid() const {return dtDid_;} ///< returns the last performed timestep

  double getDtTry(            ) const {return dtTry_;} ///< returns the timestep to try in the next step

  /// Sets the timestep to try in the next step
  /**
   * This should be used if there is some factor besides the normal ODE evolution that may modify the timestep (e.g. a quantum jump in a Monte Carlo wave-function evolution).
   * For the most part, this will be a timestep *decrease*.
   */
  void   setDtTry(double dtTry)       {dtTry_=dtTry;}

  /// \name Getters of precision parameters
  //@{
  double getEpsRel() const {return epsRel_;} ///< relative precision
  double getEpsAbs() const {return epsAbs_;} ///< absolute precision
  //@}

  void update(double t, double dtTry);

  TimeStepBookkeeper(const TimeStepBookkeeper&);
  TimeStepBookkeeper& operator=(const TimeStepBookkeeper&); ///< straightforward assignment operator that avoids self-assignment

protected:
  /// straightforward constructor
  TimeStepBookkeeper(double dtInit, ///< the initial timestep to try in the first step
                     double epsRel, ///< relative precision
                     double epsAbs  ///< absolute precision
                    );

private:
#ifdef BZ_HAVE_BOOST_SERIALIZATION
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive& ar, const unsigned int) {ar & t_ & dtDid_ & dtTry_;}
#endif // BZ_HAVE_BOOST_SERIALIZATION

  double t_, dtTry_, dtDid_;

  const double epsRel_, epsAbs_;

};


/// Template-parameter independent base class taking care of logging
class LoggingBase
{
public:
  LoggingBase() : nDerivsCalls_(), nSteps_(), nFailedSteps_() {}
  
  size_t nDerivsCalls() const {return nDerivsCalls_;}
  size_t nSteps      () const {return nSteps_      ;}
  size_t nFailedSteps() const {return nFailedSteps_;}

  void registerDerivsCall () const {++nDerivsCalls_;}
  void registerStep       ()       {++nSteps_      ;}
  void registerFailedSteps(size_t failedStepsLast) {nFailedSteps_+=failedStepsLast;}
  
  std::ostream& logOnEnd(std::ostream& os) const;

private:
#ifdef BZ_HAVE_BOOST_SERIALIZATION
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive& ar, const unsigned int) {ar & nDerivsCalls_ & nSteps_ & nFailedSteps_;}
#endif // BZ_HAVE_BOOST_SERIALIZATION

  mutable size_t nDerivsCalls_;
  size_t nSteps_, nFailedSteps_;

};


/// Class for serialization of Evolved states
/**
 * \see trajectory::AdaptiveIO
 *
 * \tparam A the array type
 *
 * \todo Think about using a shared pointer instead of plain reference for referencing the array
 *
 */
template<typename A>
class EvolvedIO : public TimeStepBookkeeper, public LoggingBase
{
  EvolvedIO(const EvolvedIO&) = delete; EvolvedIO& operator=(const EvolvedIO&) = delete;

public:
  typedef std::shared_ptr<EvolvedIO>            Ptr;
  typedef std::shared_ptr<const EvolvedIO> ConstPtr;

  /// straightforward constructor \see TimeStepBookkeeper::TimeStepBookkeeper()
  EvolvedIO(A&& a, double dtInit, double epsRel, double epsAbs)
    : TimeStepBookkeeper{dtInit,epsRel,epsAbs}, LoggingBase(), a_{std::forward<A>(a)}
  {}

  using TimeStepBookkeeper::operator=;

  A      & getA()       {return this->a_;}
  A const& getA() const {return this->a_;}

  virtual ~EvolvedIO() {} ///< necessary in order that EvolvedIO be polymorphic

private:
#ifdef BZ_HAVE_BOOST_SERIALIZATION

  /// The serialization of A by reference leads to memory leak (for not completely understood reasons),
  /** hence we adopt serialization by a temporary, which necessitates splitting save/load. */
  friend class boost::serialization::access;
  template<class Archive>
  void save(Archive& ar, const unsigned int) const
  {A temp(a_); ar & temp & boost::serialization::base_object<LoggingBase>(*this) & boost::serialization::base_object<TimeStepBookkeeper>(*this);}

  template<class Archive>
  void load(Archive& ar, const unsigned int)
  {A temp; ar & temp & boost::serialization::base_object<LoggingBase>(*this) & boost::serialization::base_object<TimeStepBookkeeper>(*this); a_.reference(temp);}

  BOOST_SERIALIZATION_SPLIT_MEMBER()

#endif // BZ_HAVE_BOOST_SERIALIZATION

  A a_;

};



/// A common interface for (adaptive stepsize) ODE drivers
/**
 * It takes the array type it operates on as template parameter. A given array type can be adapted to the form expected by Evolved by suitable specializations of “memory traits” functions.
 * (cf. ArrayTraits.h, for an implementation for `blitz::Array` cf. BlitzArray.h)
 *
 * The array which is actually "evolved" is not owned by Evolved.
 *
 * The class uses the strategy idiom for calculating the time derivative. The use of [<tt>boost::function</tt>](http://www.boost.org/doc/libs/1_55_0/doc/html/function.html)
 * assures that a rather wide range of entities will be accepted as strategy functor.
 *
 * \tparam A the array type
 *
 * \see MakerGSL::_
 *
 */
template<typename A>
class Evolved : public EvolvedIO<A>
{
public:
  typedef std::function<void(double, const A&, A&)> Derivs; ///< the strategy functor to calculate time derivative at a given time (3rd argument for output)

  typedef std::shared_ptr<      Evolved>      Ptr;
  typedef std::shared_ptr<const Evolved> ConstPtr;
  
  /// straightforward constructor \see EvolvedIO::EvolvedIO()
  Evolved(A&& a, Derivs derivs, double dtInit, double epsRel, double epsAbs) 
    : EvolvedIO<A>(std::forward<A>(a),dtInit,epsRel,epsAbs), derivs_(derivs), countedDerivs_([this](double t, const A& y, A& dydt) {this->countedDerivs(t,y,dydt);})
  {}

  using TimeStepBookkeeper::operator=;
  using EvolvedIO<A>::getA;

  virtual ~Evolved() {}

  /// takes a single adaptive step
  void step(double deltaT, ///< *maximum* length of the timestep
            std::ostream& logStream)
  {
    if (mathutils::sign(deltaT)!=mathutils::sign(EvolvedIO<A>::getDtTry())) {
      // Stepping backward
      this->setDtTry(-EvolvedIO<A>::getDtDid());
      step_v(deltaT,logStream);
    }
    else step_v(deltaT,logStream);
    this->registerStep();
    this->registerFailedSteps(nFailedStepsLast_v());
  }
  
  std::ostream& streamParameters(std::ostream& os) const {return streamParameters_v(os);} ///< delegates to private virtual

  /// Wrapper around the derivatives that registers the number of calls
  void countedDerivs(double t, const A& y, A& dydt) const {derivs_(t,y,dydt); this->registerDerivsCall();}

  /// \name Getter
  //@{
  const Derivs getDerivs() const {return countedDerivs_;}
  //@}

  size_t nFailedStepsLast() const {return nFailedStepsLast_v();} ///< number of failed steps in the last timestep (delegates to pure virtual)
  
private:
  virtual void step_v(double deltaT, std::ostream&) = 0;
  virtual std::ostream& streamParameters_v(std::ostream&) const = 0;
  virtual size_t nFailedStepsLast_v() const = 0;

  const Derivs derivs_, countedDerivs_;

};


/// \name Generic evolution functions
//@{
/// evolves for exactly time `deltaT`
/** \tparam E type of the object to evolve. Implicit interface assumed: member function named step with signature `...(double)` */
template<typename E>
void evolve(E& e, double deltaT, std::ostream& logStream=std::clog)
{
  double endTime=e.getTime()+deltaT;
  while (double dt=endTime-e.getTime()) e.step(dt,logStream);
}


/// evolves up to exactly time `t` \copydetails evolve
template<typename E>
void evolveTo(E& e, double t, std::ostream& logStream=std::clog)
{
  evolve(e,t-e.getTime(),logStream);
}
//@}



/// Factory class for Evolved types
/** \tparam A the array type for Evolved */
template<typename A>
class Maker
{
public:
  typedef typename Evolved<A>::Ptr Ptr;
  typedef typename Evolved<A>::Derivs Derivs;
  
  /// The factory member function expecting the most generic set of parameters
  template <typename ARRAY>
  const Ptr operator()(ARRAY&& array, Derivs derivs, double dtInit, double epsRel, double epsAbs,
                       const A& scaleAbs ///< this parameter is basically an element-wise mask for stepsize control – for an explanation cf. [<tt>gsl_odeiv_control_scaled_new</tt>](http://www.gnu.org/software/gsl/manual/html_node/Adaptive-Step_002dsize-Control.html#Adaptive-Step_002dsize-Control)
                      ) const {return make(std::move(array),derivs,dtInit,epsRel,epsAbs,scaleAbs);}

  virtual ~Maker() {}
  
private:
  virtual const Ptr make(A&&, Derivs, double dtInit, double epsRel, double epsAbs, const A& scaleAbs) const = 0;

};


template<typename A>
auto makeIO(A&& a, double time=0)
{
  typename EvolvedIO<A>::Ptr res = std::make_shared<EvolvedIO<A>>(std::forward<A>(a),0,0,0);
  res->setTime(time);
  return res;
}


} // evolved


#endif // CPPQEDCORE_UTILS_EVOLVED_H_INCLUDED
