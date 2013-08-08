// -*- C++ -*-
#ifndef STRUCTURE_AVERAGED_H_INCLUDED
#define STRUCTURE_AVERAGED_H_INCLUDED

#include "AveragedFwd.h"

#include "LiouvilleanAveragedCommon.h"

#include "Exception.h"

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
  typedef boost::shared_ptr<const AveragedCommon> Ptr;

  typedef LiouvilleanAveragedCommon::DArray1D Averages; ///< The 1D real array storing the calculated quantum averages (perhaps in real-imaginary pairs if a given average is complex).

  virtual ~AveragedCommon() {}

  /// This function is a hook between LiouvilleanAveragedCommonRanked::average and display.
  /**
   * It allows for arbitrary manipulations with the already calculated set of averages before displaying them.
   * 
   * The most common usage is for calculation of higher moments of a quantum operator \f$\Ob\f$. In this case, since the Averages array returned by
   * LiouvilleanAveragedCommonRanked::average must depend linearly on \f$\rho,\f$ it may contain \f$\avr{\Ob}=\Tr{\Ob\rho}\f$ and \f$\avr{\Ob^2}=\Tr{\Ob^2\rho},\f$
   * but it cannot contain e.g. the variance \f$\avr{\Ob^2}-\avr{\Ob}^2\f$ because this depends quadratically on \f$\rho\f$. The variance can be calculated by this function.
   * 
   */
  void process(Averages& averages) const {process_v(averages);}

  /// Displays the system characteristics in a nicely tabulated format
  /**
   * \link LiouvilleanAveragedCommonRanked::average Calculated\endlink (& \link process processed\endlink) quantum averages (or higher moments or any other system characteristics
   * that can be calculated on the basis of quantum averages) are displayed.
   * 
   * The display format greatly depends on what kind of system is in question (in particular, \link ElementAveraged elementary system\endlink or Composite).
   * 
   */
  std::ostream& display(const Averages& averages, std::ostream& os, int i) const {display_v(averages,os,i); return os;}

private:
  virtual void process_v(      Averages&                    ) const = 0;
  virtual void display_v(const Averages&, std::ostream&, int) const = 0;

};


///////////
//
// Averaged
//
///////////


/// The first partial specialization of the general template Averaged for the one-time dependence case (\link TimeDependence Cases 1 & 2\endlink)
/** It simply composes LiouvilleanAveragedCommonRanked and AveragedCommon. */
template<int RANK>
class Averaged<RANK,true>
  : public quantumdata::Types<RANK,LiouvilleanAveragedCommonRanked<RANK> >,public AveragedCommon
{
public:
  static const int N_RANK=RANK;

  typedef boost::shared_ptr<const Averaged> Ptr;

  typedef AveragedCommon::Averages Averages;

  typedef quantumdata::LazyDensityOperator<RANK> LazyDensityOperator;

};


/// The second partial specialization of the general template Averaged for the no-time dependence case (\link TimeDependence Cases 3 & 4\endlink)
/** Similarly to Liouvillean<RANK,false>, it simply forwards the time-dependent virtual functions to time-independent ones. */
template<int RANK>
class Averaged<RANK,false>
  : public Averaged<RANK,true>
{
public:
  typedef typename Averaged<RANK,true>::Averages            Averages           ;
  typedef typename Averaged<RANK,true>::LazyDensityOperator LazyDensityOperator;

private:
  const Averages average_v(double, const LazyDensityOperator& matrix) const {return average_v(matrix);}

  virtual const Averages average_v(const LazyDensityOperator&) const = 0;

};


} // structure

#endif // STRUCTURE_AVERAGED_H_INCLUDED
