// -*- C++ -*-
#ifndef   _MASTER_FWD_H
#define   _MASTER_FWD_H

#include "TMP_Tools.h"


namespace quantumtrajectory {

namespace master {

struct NonUnitaryIP ;
struct NoLiouvillean;

template<int RANK>
class Base;

} // master


template<int RANK, typename V=tmptools::V0, bool IS_FAST=false>
class Master;

} // quantumtrajectory


#endif // _MASTER_FWD_H
