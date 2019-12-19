// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef   CPPQEDCORE_QUANTUMTRAJECTORY_PROJECTINGMCWF_TRAJECTORY_H_INCLUDED
#define   CPPQEDCORE_QUANTUMTRAJECTORY_PROJECTINGMCWF_TRAJECTORY_H_INCLUDED

#include "ProjectingMCWF_TrajectoryFwd.h"

#include "MCWF_Trajectory.h"

#include <boost/ptr_container/ptr_vector.hpp>


namespace quantumtrajectory {


/// Derived from MCWF_Trajectory, this class uses a set of reference state-vectors to project the evolved state vector on
/**
 * The aim is to compare the evolved state against a certain set (eventually, a complete basis) of reference states.
 * This might form a nonorthogonal set, so that we adopt the \ref nonorthogonalformalism covariant-contravariant formalism.
 * 
 * \see This possibility was first developed for our paper Ref. \cite vukics2009cavity
 * 
 * The projections (the elements of the state vector in the given basis – that might not add up to 1 if the basis is incomplete) are displayed
 * in addition to each \link MCWF_Trajectory::display_v trajectory display\endlink. The \link MCWF_Trajectory::displayKey_v key\endlink is supplemented accordingly.
 * 
 */
template<int RANK>
class ProjectingMCWF_Trajectory : public MCWF_Trajectory<RANK>
{
private:
  typedef MCWF_Trajectory<RANK> Base;

  typedef typename Base::StateVector    StateVector   ;
  typedef typename Base::StateVectorLow StateVectorLow;

public:
  typedef boost::ptr_vector<StateVector> Basis; ///< The set of reference states to compare againts is stored as a \refBoostConstruct{ptr_vector,ptr_container/doc/ptr_vector.html}

  /// The signature is identical to MCWF_Trajectory::MCWF_Trajectory, but the Basis set must be supplied as well.
  template<typename SYS>
  ProjectingMCWF_Trajectory(
                            StateVector& psi,
                            const Basis& basis, ///< the set of reference states to compare the MCWF-evolved state vector against
                            const SYS& sys,
                            const mcwf::Pars& p,
                            const StateVectorLow& scaleAbs=StateVectorLow()
                            )
    : Base(psi,sys,p,scaleAbs), basis_(basis), metricTensor_uu_(help())
  {}

private:
  std::ostream&    display_v(std::ostream&, int    ) const override;
  std::ostream& displayKey_v(std::ostream&, size_t&) const override;

  const linalg::CMatrix help() const;

  const Basis basis_;
  const linalg::CMatrix metricTensor_uu_;

};


} // quantumtrajectory



#endif // CPPQEDCORE_QUANTUMTRAJECTORY_PROJECTINGMCWF_TRAJECTORY_H_INCLUDED
