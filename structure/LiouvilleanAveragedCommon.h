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
  static const DArray1D defaultArray;

  virtual ~LiouvilleanAveragedCommon() {}

  // For a  Liouvillean, the key is a description of decay channels, for an Averaged, it describes the displayed columns
  static void displayKey(std::ostream& o, size_t& i, Ptr common)
  {
    if (common) common->displayKey(o,i);
  }

  virtual void displayKey(std::ostream&, size_t&) const = 0;

};


} // structure



#endif // STRUCTURE_LIOUVILLEANAVERAGEDCOMMON_H_INCLUDED
