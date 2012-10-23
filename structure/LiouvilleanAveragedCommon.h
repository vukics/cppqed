// -*- C++ -*-
#ifndef   STRUCTURE_LIOUVILLEANAVERAGEDCOMMON_H_INCLUDED
#define   STRUCTURE_LIOUVILLEANAVERAGEDCOMMON_H_INCLUDED

#include "LiouvilleanAveragedCommonFwd.h"

#include "BlitzArray.h"

#include <boost/shared_ptr.hpp>

namespace structure {


class LiouvilleanAveragedCommon
{
public:
  typedef boost::shared_ptr<const LiouvilleanAveragedCommon> Ptr;

  typedef TTD_DARRAY(1) DArray1D;

  virtual ~LiouvilleanAveragedCommon() {}

  // For a  Liouvillean, the key is a description of decay channels, for an Averaged, it describes the displayed columns
  void displayKey(std::ostream& o, size_t& i) const {displayKey_v(o,i);}

private:
  virtual void displayKey_v(std::ostream&, size_t&) const = 0;

};


} // structure



#endif // STRUCTURE_LIOUVILLEANAVERAGEDCOMMON_H_INCLUDED
