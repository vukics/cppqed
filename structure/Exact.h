// -*- C++ -*-
/// \briefFile{Defines the hierarchical partial specializations of structure::Exact}
#ifndef STRUCTURE_EXACT_H_INCLUDED
#define STRUCTURE_EXACT_H_INCLUDED

#include "ExactFwd.h"

#include "Types.h"

#include <boost/shared_ptr.hpp>

namespace structure {


/// The template-parameter-independent base of Exact
class ExactCommon
{
public:
  typedef boost::shared_ptr<const ExactCommon> Ptr;

  virtual ~ExactCommon() {}

  bool isUnitary() const {return isUnitary_v();} ///< Describes whether the interaction picture is unitary

private:
  virtual bool isUnitary_v() const = 0;

};


/// The first partial specialization of the general template Exact for the most general, two-time dependence 
/** This corresponds to \link TimeDependenceLevel Case 1\endlink (#TWO_TIME).*/
template<int RANK>
class Exact<RANK,true> : public ExactCommon, private quantumdata::Types<RANK> 
{
public:
  typedef boost::shared_ptr<const Exact> Ptr;

  typedef typename quantumdata::Types<RANK>::StateVectorLow StateVectorLow;

  virtual ~Exact() {}

  /// Describes the operation which transforms from interaction picture to the normal picture: \f$\ket{\Psi(t)}\rightarrow U(t,t_0)\ket{\Psi}\f$
  void actWithU(double t, ///<[in] \f$t\f$
                StateVectorLow& psi, ///<[in/out] \f$\ket{\Psi}\f$
                double tIntPic0 ///<[in] \f$t_0\f$
               ) const {return actWithU_v(t,psi,tIntPic0);} 

private:
  virtual void actWithU_v(double, StateVectorLow&, double) const = 0;

};


/// The second partial specialization of the general template Exact for one-time dependence 
/** This corresponds to \link TimeDependenceLevel Case 3\endlink of the two possibilities for #ONE_TIME. */
template<int RANK>
class Exact<RANK,false> : public Exact<RANK,true>
{
public:
  typedef typename Exact<RANK,true>::StateVectorLow StateVectorLow;
  
private:
  void actWithU_v(double t, StateVectorLow& psi, double tIntPic0) const {actWithU_v(t-tIntPic0,psi);}

  virtual void actWithU_v(double t, StateVectorLow& psi) const = 0;

};


} // structure

#endif // STRUCTURE_EXACT_H_INCLUDED
