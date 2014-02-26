// -*- C++ -*-
#ifndef   CPPQEDCORE_UTILS_STOCHASTICTRAJECTORYFWD_H_INCLUDED
#define   CPPQEDCORE_UTILS_STOCHASTICTRAJECTORYFWD_H_INCLUDED


namespace trajectory {

struct ParsStochastic;

template<typename T> 
class Averageable;

template<typename T, typename T_ELEM=T>
class Ensemble;

namespace ensemble {

template<typename T, typename T_ELEM>
class Traits;

}

template<typename A, typename T> 
class Stochastic;


}

#endif // CPPQEDCORE_UTILS_STOCHASTICTRAJECTORYFWD_H_INCLUDED
