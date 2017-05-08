// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
/// \briefFile{Defines tools related to the description of different time-dependence levels}
#ifndef   CPPQEDCORE_STRUCTURE_TIME_H_INCLUDED
#define   CPPQEDCORE_STRUCTURE_TIME_H_INCLUDED

#include "TimeFwd.h"

#include <boost/operators.hpp>

#include <boost/mpl/if.hpp>


namespace structure {


/// Enumeration of different possibilities for time dependence of Hamiltonians
/** With \f$t_0\f$ being the time instant where the \link Exact two pictures\endlink coincide (cf. quantumtrajectory::QuantumTrajectory::getT0): */
enum TimeDependence {
  TWO_TIME, ///< **Case 1** \f$H(t,t_0)\f$ – Time-dependent problem + \link Exact exact part\endlink (\f$U(t,t_0)\f$)
  ONE_TIME, ///< **Case 2** \f$H(t)\f$ – Time-dependent problem, no exact part **OR Case 3** \f$H(t-t_0)\f$ – Time-independent problem + \link Exact exact part\endlink (\f$U(t-t_0)\f$)
  NO_TIME   ///< **Case 4** \f$H(0)\f$ – Time-independent problem, no exact part
};


/// Describes two-time dependence corresponding to the case #TWO_TIME
class TwoTime : private boost::equality_comparable<TwoTime>
{
public:
  TwoTime(double t, double t0) : t_(t), t0_(t0) {}
  
  double getT () const {return t_ ;}
  double getT0() const {return t0_;}
  
private:
  friend bool operator==(TwoTime t1, TwoTime t2) {return t1.t_==t2.t_ && t1.t0_==t2.t0_;}
  
  double t_, t0_;
  
};

/// Describes one-time dependence corresponding to the case #ONE_TIME
class OneTime : private boost::equality_comparable<OneTime>
{
public:
  explicit OneTime(double t, double t0=0) : deltaT_(t-t0) {}

  operator double() const {return deltaT_;}

private:
  friend bool operator==(OneTime t1, OneTime t2) {return t1.deltaT_==t2.deltaT_;}

  double deltaT_;
  
};

/// Describes no-time dependence corresponding to the case #NO_TIME
class NoTime
{
public:
  explicit NoTime(double, double=0) {}

};


/// Comprises tools related to the description of different time-dependence levels
/** \see HamiltonianTimeDependenceDispatched for a paradigmatic example of usage of the time bundle */
namespace time {


/// Metafunction dispatching the three classes TwoTime, OneTime, & NoTime according to the template parameter `TD`
template<TimeDependence TD>
struct Dispatcher : boost::mpl::if_c<TD==TWO_TIME,TwoTime,
                                     typename boost::mpl::if_c<TD==ONE_TIME,OneTime,NoTime>::type
                                     >
{};


/// Metafunction dispatching two OneTime & NoTime according to the template parameter `IS_TIME_DEPENDENT`
/** \see LiouvilleanTimeDependenceDispatched for a paradigmatic example of usage */
template<bool IS_TIME_DEPENDENT>
struct DispatcherIsTimeDependent : boost::mpl::if_c<IS_TIME_DEPENDENT,OneTime,NoTime>
{};


/// Metafunction dispatching two TwoTime & OneTime according to the template parameter `IS_TWO_TIME`
/** \see ExactTimeDependenceDispatched for a paradigmatic example of usage */
template<bool IS_TWO_TIME>
struct DispatcherIsTwoTime : boost::mpl::if_c<IS_TWO_TIME,TwoTime,OneTime>
{};

} // time


} // structure

#endif // CPPQEDCORE_STRUCTURE_TIME_H_INCLUDED
