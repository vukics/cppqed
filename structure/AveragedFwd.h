// -*- C++ -*-
#ifndef   STRUCTURE_AVERAGEDFWD_H_INCLUDED
#define   STRUCTURE_AVERAGEDFWD_H_INCLUDED

namespace structure {

class AveragedCommon;


/// The interface every system that calculates and displays quantum averages must present towards the trajectory drivers
/**
 * \tparamRANK
 * \tparam IS_TIME_DEPENDENT describes whether the operators whose quantum average is calculated are time-dependent. Default `true`, the most general case.
 * 
 *  Similarly to Exact & Liouvillean, the general template is never defined, but a hierarchy of partial specializations in the second template argument.
 * 
 */
template<int, bool IS_TIME_DEPENDENT=true> class Averaged;

} // structure


#endif // STRUCTURE_AVERAGEDFWD_H_INCLUDED
