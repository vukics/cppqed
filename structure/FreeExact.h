// -*- C++ -*-
/// \briefFileDefault
#ifndef STRUCTURE_FREEEXACT_H_INCLUDED
#define STRUCTURE_FREEEXACT_H_INCLUDED

#include "Exact.h"

#include <boost/utility/enable_if.hpp>


namespace structure {

/// A unary implementation of Exact assuming that the operator that transforms between the pictures is diagonal
/**
 * If this is not the case, then the system has to inherit directly from Exact<1> and implement the Exact::actWithU function.
 * 
 * \tparam IS_TWO_TIME Same meaning as for Exact, but here the default is `false`
 * 
 */
template<bool IS_TWO_TIME=false>
class FreeExact;



namespace details {


class FreeExactBase
{
public:
  typedef CArray<1> Diagonal;
  
protected:
  explicit FreeExactBase(size_t dim) : diagonal_(int(dim)) {}

  Diagonal& getDiagonal() const {return diagonal_;}
  
private:
  mutable Diagonal diagonal_;

};


} // details


/// Specialization of FreeExact for two-time dependence case
template<>
class FreeExact<true> : public Exact<1,true>, public details::FreeExactBase
{
protected:
  explicit FreeExact(size_t dim) : details::FreeExactBase(dim), t_(0), t0_(0) {}
  
private:
  void actWithU_v(double t, StateVectorLow& psi, double t0) const {if (t!=t_ || t0!=t0_) {updateU(t_=t,t0_=t0);} psi*=getDiagonal();} ///< Implements Exact<1,true >
  
  virtual void updateU(double,double) const = 0; ///< Updates diagonals to the given \f$t\f$ & \f$t_0\f$
  
private:
  mutable double t_, t0_;
};


/// Specialization of FreeExact for one-time dependence case
template<>
class FreeExact<false> : public Exact<1,false>, public details::FreeExactBase
{
protected:
  explicit FreeExact(size_t dim) : details::FreeExactBase(dim), deltaT_(0) {}

private:
  void actWithU_v(double deltaT, StateVectorLow& psi) const {if (deltaT!=deltaT_) {updateU(deltaT_=deltaT);} psi*=getDiagonal();} ///< Implements Exact<1,false>

  virtual void updateU(double deltaT) const = 0; ///< Updates diagonals to the given \f$t-t_0\f$

  mutable double deltaT_;

};

} // structure


#endif // STRUCTURE_FREEEXACT_H_INCLUDED
