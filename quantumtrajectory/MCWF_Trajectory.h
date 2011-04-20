// -*- C++ -*-
#ifndef _MCWF_TRAJECTORY_H
#define _MCWF_TRAJECTORY_H

#include "MCWF_TrajectoryFwd.h"

#include "DimensionsBookkeeperFwd.h"
#include "QuantumSystemFwd.h"
#include "ExactFwd.h"
#include "HamiltonianFwd.h"
#include "LiouvilleanFwd.h"
#include "AveragedFwd.h"

#include "Types.h"

#include "StochasticTrajectory.h"

#include<boost/tuple/tuple.hpp>

#include<boost/mpl/identity.hpp>


namespace mpl=boost::mpl;

namespace quantumtrajectory {


class MCWF_TrajectoryFileOpeningException : public cpputils::TaggedException
{
public:
  MCWF_TrajectoryFileOpeningException(const std::string tag) : cpputils::TaggedException(tag) {}

};

///////////////////////////////////////
//
// Monte Carlo wave-function trajectory
//
///////////////////////////////////////

// NEEDS_WORK factor out template-parameter independent code




#define BASE_class trajectory::StochasticTrajectory<typename quantumdata::Types<RANK>::StateVectorLow, const quantumdata::StateVector<RANK>&>

template<int RANK>
class MCWF_Trajectory : public BASE_class
{
public:
  typedef quantumdata::StateVector<RANK> StateVector;

  typedef typename StateVector::StateVectorLow StateVectorLow;

  typedef structure::QuantumSystem<RANK> QuantumSystem;
  typedef structure::Exact        <RANK> Exact        ;
  typedef structure::Hamiltonian  <RANK> Hamiltonian  ;
  typedef structure::Liouvillean  <RANK> Liouvillean  ;
  typedef structure::Averaged     <RANK> Averaged     ;  

  typedef BASE_class Base;

#undef  BASE_class

  typedef boost::tuple<int,StateVectorLow> indexSVL_tuple;

  using Base::getEvolved; using Base::getRandomized; using Base::getOstream; using Base::getDtDid; using Base::getTime;


  MCWF_Trajectory(
		  StateVector& psi,
		  const QuantumSystem& sys,
		  const ParsMCWF_Trajectory&,
		  const StateVectorLow& =StateVectorLow()
		  );

  virtual ~MCWF_Trajectory();

  void derivs(double, const StateVectorLow&, StateVectorLow&) const;

  void step(double) const; // performs one single adaptive-stepsize MCWF step of specified maximal length

  void   displayParameters(   ) const;

protected:
  virtual size_t displayMoreKey (   ) const;

  virtual void   displayEvenMore(int) const {}

  const StateVector& getPsi      () const {return psi_;} 
  const StateVector& toBeAveraged() const {return psi_;} 

private:
  void displayMore(int) const;

  mutable double tIntPic0_ ; // The time instant of the beginning of the current time step.
  // mutable double random_   ; 
  // Maybe put into StochasticTrajectory somehow.

  StateVector& psi_;

  const QuantumSystem*const qs_;
  const Exact        *const ex_; 
  const Hamiltonian  *const ha_;
  const Liouvillean  *const li_; 
  const Averaged     *const av_;

  const double dpLimit_;

  const unsigned svdc_;

  const std::string file_;

  const std::string initFile_;

};


} // quantumtrajectory


#include "impl/MCWF_Trajectory.tcc"


#endif // _MCWF_TRAJECTORY_H
