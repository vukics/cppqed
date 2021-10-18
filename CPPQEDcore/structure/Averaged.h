// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_STRUCTURE_AVERAGED_H_INCLUDED
#define CPPQEDCORE_STRUCTURE_AVERAGED_H_INCLUDED

#include "LiouvilleanAveragedCommon.h"
#include "TimeDependence.h"

#include <iosfwd>


namespace structure {


/////////////////
//
// AveragedCommon
//
/////////////////


/// The template-parameter independent base of Averaged
class AveragedCommon
{
public:
  virtual ~AveragedCommon() {}

  /// This function is a hook between LiouvilleanAveragedCommonRanked::average and stream.
  /**
   * It allows for arbitrary manipulations with the already calculated set of averages before streaming them.
   * 
   * The most common usage is for calculation of higher moments of a quantum operator \f$\Ob\f$. In this case, since the Averages array returned by
   * LiouvilleanAveragedCommonRanked::average must depend linearly on \f$\rho,\f$ it may contain \f$\avr{\Ob}=\Tr{\Ob\rho}\f$ and \f$\avr{\Ob^2}=\Tr{\Ob^2\rho},\f$
   * but it cannot contain e.g. the variance \f$\avr{\Ob^2}-\avr{\Ob}^2\f$ because this depends quadratically on \f$\rho\f$. The variance can be calculated by this function.
   * 
   */
  void process(Averages& averages) const {process_v(averages);}

  /// Streams the system characteristics in a nicely tabulated format
  /**
   * \link LiouvilleanAveragedCommonRanked::average Calculated\endlink (& \link process processed\endlink) quantum averages (or higher moments or any other system characteristics
   * that can be calculated on the basis of quantum averages) are streamed.
   * 
   * The display format greatly depends on what kind of system is in question (in particular, \link ElementAveraged elementary system\endlink or Composite).
   * 
   * \todo The precision argument is obsolete after rev. #04265ae since precision might be propagated just via FormDouble::overallPrecision throughout the whole framework.
   *
   */
  std::ostream& stream(const Averages& averages, std::ostream& os, int precision) const {return stream_v(averages,os,precision);}

private:
  virtual void          process_v(      Averages&                    ) const {}
  virtual std::ostream& stream_v(const Averages&, std::ostream&, int) const = 0;

};


///////////
//
// Averaged
//
///////////


/// The interface every system that calculates and streams quantum averages must present towards the trajectory drivers
/**
 * It simply composes LiouvilleanAveragedCommonRanked and AveragedCommon.
 * 
 * \tparamRANK
 * 
 * \see The design is exactly the same as that of Liouvillean
 * 
 */
template<int RANK>
class Averaged : public LiouvilleanAveragedCommonRanked<RANK>, public AveragedCommon
{
public:
  static const int N_RANK=RANK;

};


template <int RANK>
using AveragedPtr=std::shared_ptr<const Averaged<RANK>>;


/// Implements the general Liouvillean interface by dispatching the two possible \link time::DispatcherIsTimeDependent time-dependence levels\endlink
/**
 * Similarly to LiouvilleanTimeDependenceDispatched, it simply forwards the time-dependent virtual functions to time-independence-dispatched ones.
 * 
 * \tparamRANK
 * \tparam IS_TIME_DEPENDENT describes whether the observables whose quantum average we want to calculate are time-dependent. `true`: OneTime – `false`: NoTime
 * 
 */
template<int RANK, bool IS_TIME_DEPENDENT>
class AveragedTimeDependenceDispatched : public Averaged<RANK>
{
public:
  typedef time::DispatcherIsTimeDependent_t<IS_TIME_DEPENDENT> Time;

private:
  const Averages average_v(double t, const quantumdata::LazyDensityOperator<RANK>& matrix) const final {return average_v(Time(t),matrix);}

  virtual const Averages average_v(Time, const quantumdata::LazyDensityOperator<RANK>&) const = 0;

};


} // structure

#endif // CPPQEDCORE_STRUCTURE_AVERAGED_H_INCLUDED
