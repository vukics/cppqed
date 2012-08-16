// -*- C++ -*-
#ifndef   UTILS_INCLUDE_STOCHASTICTRAJECTORYFWD_H_INCLUDED
#define   UTILS_INCLUDE_STOCHASTICTRAJECTORYFWD_H_INCLUDED


namespace trajectory {

struct ParsStochasticTrajectory;

template<typename T> 
class StochasticTrajectoryBase;

template<typename T, typename T_ELEM=T>
class EnsembleTrajectories;

template<typename T, typename T_ELEM>
class EnsembleTrajectoriesTraits;


template<typename A, typename T> 
class StochasticTrajectory;


}

#endif // UTILS_INCLUDE_STOCHASTICTRAJECTORYFWD_H_INCLUDED
