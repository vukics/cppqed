// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
/// \briefFileDefault
#ifndef CPPQEDCORE_STRUCTURE_FREEEXACT_H_INCLUDED
#define CPPQEDCORE_STRUCTURE_FREEEXACT_H_INCLUDED

#include "Exact.h"

#include <boost/utility/enable_if.hpp>


namespace structure {

/// A unary implementation of Exact assuming that the operator that transforms between the pictures is diagonal
/**
 * If this is not the case, then the system has to inherit directly from Exact<1> and implement the Exact::actWithU function.
 * 
 * \tparam IS_TWO_TIME Same meaning as for ExactTimeDependenceDispatched
 * 
 */
template<bool IS_TWO_TIME>
class FreeExact : public ExactTimeDependenceDispatched<1,IS_TWO_TIME>
{
public:
  typedef CArray<1> Diagonal;

  typedef typename ExactTimeDependenceDispatched<1,IS_TWO_TIME>::StateVectorLow StateVectorLow;
  
  typedef typename ExactTimeDependenceDispatched<1,IS_TWO_TIME>::Time Time;

protected:
  explicit FreeExact(size_t dim) : diagonal_(int(dim)), t_(0,0) {}
  
  Diagonal& getDiagonal() const {return diagonal_;}

private:
  virtual void actWithU_v(Time t, StateVectorLow& psi) const final {if (t!=t_) {updateU(t_=t);} psi*=getDiagonal();} ///< Implements Exact<1,true >
  
  virtual void updateU(Time) const = 0; ///< Updates diagonals to the given \f$t\f$ & \f$t_0\f$
  
private:
  mutable Diagonal diagonal_;
  mutable Time t_;

};


} // structure


#endif // CPPQEDCORE_STRUCTURE_FREEEXACT_H_INCLUDED
