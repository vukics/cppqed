// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef   CPPQEDCORE_STRUCTURE_LIOUVILLEANAVERAGEDCOMMON_H_INCLUDED
#define   CPPQEDCORE_STRUCTURE_LIOUVILLEANAVERAGEDCOMMON_H_INCLUDED

#include "LiouvilleanAveragedCommonFwd.h"

#include "Types.h"

#include "BlitzArray.h"
#include "BlitzArrayExtensions.h"
#include "Exception.h"

#include <boost/shared_ptr.hpp>

namespace structure {

/// Enumeration describing a (compile-time) choice between Liouvillean & Averaged
enum LiouvilleanAveragedTag {
  LA_Li, ///< Liouvillean
  LA_Av  ///< Averaged
};

/// A tagging struct corresponding to LiouvilleanAveragedTag
template<LiouvilleanAveragedTag>
struct LiouvilleanAveragedTag_ {LiouvilleanAveragedTag_() {}};

typedef LiouvilleanAveragedTag_<LA_Li> LA_Li_tagType;
typedef LiouvilleanAveragedTag_<LA_Av> LA_Av_tagType;

const LA_Li_tagType LA_Li_tag;
const LA_Av_tagType LA_Av_tag;


/// The template-parameter independent part of LiouvilleanAveragedCommonRanked
class LiouvilleanAveragedCommon
{
public:
  typedef boost::shared_ptr<const LiouvilleanAveragedCommon> Ptr;

  typedef DArray<1> DArray1D; ///< A 1D real array storing the quantum averages – even if they are complex, their real & imaginary parts must be stored separately as reals

  virtual ~LiouvilleanAveragedCommon() {}

  /// Displays a key (a.k.a. legend)
  /**
   * - for a  Liouvillean, the key is a description of decay channels,
   * - for an Averaged, it describes the displayed columns
   */
  std::ostream& displayKey(std::ostream& o,
                           size_t& i ///< the ordinal number where the key of the present element begins
                          ) const {return displayKey_v(o,i);}

  /// Returns the number of calculated quantum averages
  /**
   * - for a Liouvillean the quantum averages are the jump rates
   * - for an Averaged, they are the quantum averages calculated for display
   */
  size_t nAvr() const {return nAvr_v();}
  
private:
  virtual std::ostream& displayKey_v(std::ostream&, size_t&) const = 0;
  virtual size_t        nAvr_v      (                      ) const = 0;


};


/// Exception for LiouvilleanAveragedCommonRanked::average
struct AveragesNumberMismatchException : cpputils::Exception
{
  AveragesNumberMismatchException(int size, size_t nAvr) {std::cerr<<size<<' '<<nAvr<<std::endl;}
  
};


/// Exception for LiouvilleanAveragedCommonRanked::average
struct InfiniteDetectedException : cpputils::Exception {};


/// Common functionality of Liouvillean & Averaged
/** \tparamRANK */
template<int RANK>
class LiouvilleanAveragedCommonRanked : public LiouvilleanAveragedCommon
{
public:
  typedef boost::shared_ptr<const LiouvilleanAveragedCommonRanked> Ptr;
  
  typedef typename LiouvilleanAveragedCommon::DArray1D DArray1D;

  typedef quantumdata::LazyDensityOperator<RANK> LazyDensityOperator;

  virtual ~LiouvilleanAveragedCommonRanked() {}

  /// Calculates quantum averages & checks post-conditions
  /**
   * \warning The elements of the returned array must have a linear dependence on `matrix`, in particular, they must be of the form \f$\Tr{\Ob_m\rho},\f$
   * where \f$\Ob_m\f$ is a set of quantum operators defined on the Hilbert space of the system.
   * 
   * \throw AveragesNumberMismatchException if postcondition 1. is not met
   * \throw InfiniteDetectedException if postcondition 2. is not met
   * 
   */
  const DArray1D average(double t, ///< one or more of the operators whose quantum average is calculated, might be time-dependent
                         const LazyDensityOperator& matrix /// the state of the quantum system
                        ) const
  {
    const DArray1D averages(average_v(t,matrix));
    if (size_t(averages.size())!=nAvr()) throw AveragesNumberMismatchException(averages.size(),nAvr()); /// \post 1.: number of averages must equal LiouvilleanAveragedCommon::nAvr
    if (!all(blitzplusplus::isfinite(averages))) throw InfiniteDetectedException(); /// \post 2.: all quantum averages must be finite

    return averages;
  }

private:
  virtual const DArray1D average_v(double, const LazyDensityOperator&) const = 0;
  
};


} // structure



#endif // CPPQEDCORE_STRUCTURE_LIOUVILLEANAVERAGEDCOMMON_H_INCLUDED
