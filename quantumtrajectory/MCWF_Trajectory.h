// -*- C++ -*-
#ifndef QUANTUMTRAJECTORY_MCWF_TRAJECTORY_H_INCLUDED
#define QUANTUMTRAJECTORY_MCWF_TRAJECTORY_H_INCLUDED

#include "MCWF_TrajectoryFwd.h"

#include "MCWF_TrajectoryLogger.h"
#include "StateVectorFwd.h"
#include "Structure.h"

#include "StochasticTrajectory.h"

#include <boost/tuple/tuple.hpp>

#include <boost/mpl/identity.hpp>

#include <vector>

namespace mpl=boost::mpl;

namespace quantumtrajectory {


class MCWF_TrajectoryFileOpeningException : public cpputils::TaggedException
{
public:
  MCWF_TrajectoryFileOpeningException(const std::string tag) : cpputils::TaggedException(tag) {}

};


class MCWF_TrajectoryFileParsingException : public cpputils::TaggedException
{
public:
  MCWF_TrajectoryFileParsingException(const std::string tag) : cpputils::TaggedException(tag) {}

};


///////////////////////////////////////
//
// Monte Carlo wave-function trajectory
//
///////////////////////////////////////

// NEEDS_WORK factor out template-parameter independent code


#define BASE_class trajectory::Stochastic<typename quantumdata::Types<RANK>::StateVectorLow, const quantumdata::StateVector<RANK>&>

template<int RANK>
class MCWF_Trajectory : public BASE_class
{
public:
  typedef structure::Exact        <RANK> Exact      ;
  typedef structure::Hamiltonian  <RANK> Hamiltonian;
  typedef structure::Liouvillean  <RANK> Liouvillean;
  typedef structure::Averaged     <RANK> Averaged   ;

  typedef quantumdata::StateVector<RANK> StateVector;

  typedef typename StateVector::StateVectorLow StateVectorLow;

  typedef structure::QuantumSystemWrapper<RANK,true> QuantumSystemWrapper;

  typedef BASE_class Base;

#undef  BASE_class

  typedef boost::tuple<int,StateVectorLow> IndexSVL_tuple;

  using Base::getEvolved; using Base::getRandomized; using Base::getOstream; using Base::getPrecision; using Base::getDtDid; using Base::getDtTry; using Base::getTime;

  template<typename SYS>
  MCWF_Trajectory(
		  StateVector& psi,
		  const SYS& sys,
		  const ParsMCWF_Trajectory&,
		  const StateVectorLow& =StateVectorLow()
		  );

  virtual ~MCWF_Trajectory();

  void derivs(double, const StateVectorLow&, StateVectorLow&) const;

  const StateVector& getPsi() const {return psi_;} 

protected:
  virtual size_t displayMoreKey () const;

  virtual void   displayEvenMore() const {}

private:
  typedef std::vector<IndexSVL_tuple> IndexSVL_tuples;
  typedef typename Liouvillean::Probabilities DpOverDtSet;

  void step_v(double) const; // performs one single adaptive-stepsize MCWF step of specified maximal length

  void displayParameters_v() const;

  const StateVector& toBeAveraged_v() const {return psi_;} 

  void displayMore() const;
  
  void readState(std::ifstream &, bool onlySV=false);
  void writeState(std::ofstream &) const;
  
  double                coherentTimeDevelopment    (                                double Dt) const;
  const IndexSVL_tuples calculateDpOverDtSpecialSet(      DpOverDtSet* dpOverDtSet, double  t) const;

  bool                  manageTimeStep             (const DpOverDtSet& dpOverDtSet, evolved::TimeStepBookkeeper*, bool logControl=true) const;

  void                  performJump                (const DpOverDtSet&, const IndexSVL_tuples&, double) const;
  // helpers to step---we are deliberately avoiding the normal technique of defining such helpers, because in that case the whole MCWF_Trajectory has to be passed

  mutable double tIntPic0_ ; // The time instant of the beginning of the current time step.

  StateVector& psi_;

  const QuantumSystemWrapper qs_;

  const double dpLimit_, overshootTolerance_;

  const unsigned svdc_;
  const bool firstSVDisplay_;
  const int svdPrecision_;
  mutable long svdCount_;

  const std::string svExtension_;

  const std::string file_;

  const std::string initFile_;

  const MCWF_TrajectoryLogger logger_;

};


} // quantumtrajectory


#endif // QUANTUMTRAJECTORY_MCWF_TRAJECTORY_H_INCLUDED
