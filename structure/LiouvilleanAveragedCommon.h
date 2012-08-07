// -*- C++ -*-
#ifndef   LIOUVILLEAN_AVERAGED_COMMON_INCLUDED
#define   LIOUVILLEAN_AVERAGED_COMMON_INCLUDED

#include "LiouvilleanAveragedCommonFwd.h"

#include "BlitzArray.h"


namespace structure {


class LiouvilleanAveragedCommon
{
public:
  typedef TTD_DARRAY(1) DArray1D;
  static const DArray1D defaultArray;

  virtual ~LiouvilleanAveragedCommon() {}

  // For a  Liouvillean, the key is a description of decay channels,
  // for an Averaged   , it describes the displayed columns
  static void displayKey(std::ostream& o, size_t& i, const LiouvilleanAveragedCommon* common)
  {
    if (common) common->displayKey(o,i);
  }

  virtual void displayKey(std::ostream&, size_t&) const = 0;

};


} // structure



#endif // LIOUVILLEAN_AVERAGED_COMMON_INCLUDED
