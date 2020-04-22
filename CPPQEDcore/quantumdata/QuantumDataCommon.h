#ifndef CPPQEDCORE_QUANTUMDATA_QUANTUMDATAFWD_H_INCLUDED
#define CPPQEDCORE_QUANTUMDATA_QUANTUMDATAFWD_H_INCLUDED


namespace quantumdata {

struct ByReference {}; const ByReference byReference=ByReference();

template<int>
class StateVector;

template<int>
class DensityOperator;

template<int> 
class LazyDensityOperator;


} // quantumdata

#endif // CPPQEDCORE_QUANTUMDATA_QUANTUMDATAFWD_H_INCLUDED
