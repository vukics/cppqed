// -*- C++ -*-
#ifndef   STRUCTURE_LIOUVILLEANAVERAGEDCOMMON_H_INCLUDED
#define   STRUCTURE_LIOUVILLEANAVERAGEDCOMMON_H_INCLUDED

#include "LiouvilleanAveragedCommonFwd.h"

#include "Types.h"

#include "BlitzArray.h"
#include "BlitzArrayExtensions.h"
#include "Exception.h"

#include <boost/shared_ptr.hpp>

namespace structure {


enum LiouvilleanAveragedTag {LA_Li, LA_Av};

template<LiouvilleanAveragedTag>
struct LiouvilleanAveragedTag_ {LiouvilleanAveragedTag_() {}};

typedef LiouvilleanAveragedTag_<LA_Li> LA_Li_tagType;
typedef LiouvilleanAveragedTag_<LA_Av> LA_Av_tagType;

const LA_Li_tagType LA_Li_tag;
const LA_Av_tagType LA_Av_tag;

  
class LiouvilleanAveragedCommon
{
public:
  typedef boost::shared_ptr<const LiouvilleanAveragedCommon> Ptr;

  typedef TTD_DARRAY(1) DArray1D;

  virtual ~LiouvilleanAveragedCommon() {}

  // For a  Liouvillean, the key is a description of decay channels, for an Averaged, it describes the displayed columns
  void displayKey(std::ostream& o, size_t& i) const {displayKey_v(o,i);}

  // For a Liouvillean, nAvr is the number of jump probabilities
  size_t nAvr() const {return nAvr_v();}
  
private:
  virtual void   displayKey_v(std::ostream&, size_t&) const = 0;
  virtual size_t       nAvr_v(                      ) const = 0;


};


#ifndef   NDEBUG
struct AveragesNumberMismatchException : cpputils::Exception
{
  AveragesNumberMismatchException(int size, size_t nAvr) {std::cerr<<size<<' '<<nAvr<<std::endl;}
  
};
#endif // NDEBUG

struct InfiniteDetectedException : cpputils::Exception {};


template<int RANK>
class LiouvilleanAveragedCommonRanked : public LiouvilleanAveragedCommon
{
public:
  typedef boost::shared_ptr<const LiouvilleanAveragedCommonRanked> Ptr;
  
  typedef typename LiouvilleanAveragedCommon::DArray1D DArray1D;

  typedef quantumdata::LazyDensityOperator<RANK> LazyDensityOperator;

  virtual ~LiouvilleanAveragedCommonRanked() {}

  
  using LiouvilleanAveragedCommon::nAvr;
  
  
  const DArray1D average(double t, const LazyDensityOperator& matrix) const
  {
    const DArray1D averages(average_v(t,matrix));
    if (size_t(averages.size())!=nAvr()) throw AveragesNumberMismatchException(averages.size(),nAvr());
    if (!all(blitzplusplus::isfinite(averages))) throw InfiniteDetectedException();

    return averages;
  }

private:
  virtual const DArray1D average_v(double, const LazyDensityOperator&) const = 0;
  
};


} // structure



#endif // STRUCTURE_LIOUVILLEANAVERAGEDCOMMON_H_INCLUDED
