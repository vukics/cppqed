// -*- C++ -*-
#ifndef   QUANTUMTRAJECTORY_MASTERFWD_H_INCLUDED
#define   QUANTUMTRAJECTORY_MASTERFWD_H_INCLUDED

#include "TMP_Tools.h"


namespace quantumtrajectory {

namespace master {

struct NonUnitaryIP ;
struct NoLiouvillean;

template<int RANK>
class Base;

} // master


template<int RANK, typename V=tmptools::V_Empty, bool IS_FAST=false>
class Master;

} // quantumtrajectory


#endif // QUANTUMTRAJECTORY_MASTERFWD_H_INCLUDED
