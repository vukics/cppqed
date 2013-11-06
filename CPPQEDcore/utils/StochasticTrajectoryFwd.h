// -*- C++ -*-
#ifndef   UTILS_STOCHASTICTRAJECTORYFWD_H_INCLUDED
#define   UTILS_STOCHASTICTRAJECTORYFWD_H_INCLUDED


namespace trajectory {

struct ParsStochastic;

template<typename T> 
class Averageable;

template<typename T, typename T_ELEM=T>
class Ensemble;

template<typename T, typename T_ELEM>
class EnsembleTraits;


template<typename A, typename T> 
class Stochastic;


}

#endif // UTILS_STOCHASTICTRAJECTORYFWD_H_INCLUDED
