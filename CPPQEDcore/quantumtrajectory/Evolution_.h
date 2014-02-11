/// Highest level driver functions for quantum trajectories
#ifndef QUANTUMTRAJECTORY_EVOLUTION__H_INCLUDED
#define QUANTUMTRAJECTORY_EVOLUTION__H_INCLUDED

#include "Evolution_Fwd.h"

#include "StateVectorFwd.h"
#include "MCWF_TrajectoryFwd.h"
#include "ParsMCWF_Trajectory.h"
#include "QuantumSystemFwd.h"

#include "SmartPtr.h"
#include "TMP_Tools.h"

#include <iosfwd>



using namespace quantumtrajectory;


enum EvolutionMode {/*EM_CONVERGENCE,*/ EM_SINGLE, EM_ENSEMBLE, EM_MASTER, EM_MASTER_FAST};
// covergence runs a comparison between ensemble & single

std::ostream& operator<<(std::ostream&, EvolutionMode );
std::istream& operator>>(std::istream&, EvolutionMode&);


struct ParsEvolution : public trajectory::ParsRun, public ParsMCWF {

  EvolutionMode &evol;
  bool &negativity, &timeAverage;
  double &relaxationTime;

  ParsEvolution(parameters::ParameterTable& p, const std::string& mod="");

};


template<int RANK, typename SYS>
const boost::shared_ptr<MCWF_Trajectory<RANK> > makeMCWF(quantumdata::StateVector<RANK>&, const SYS&, const ParsEvolution&);


template<typename V, int RANK>
void evolve(quantumdata::StateVector<RANK>&, typename structure::QuantumSystem<RANK>::Ptr,
            const ParsEvolution&);



template<int RANK>
inline
void evolve(quantumdata::StateVector<RANK>& psi,
            typename structure::QuantumSystem<RANK>::Ptr sys,
            const ParsEvolution& p)
{
  evolve<tmptools::V_Empty>(psi,sys,p);
}

template<typename V, int RANK>
inline
void evolve(quantumdata::StateVector<RANK>& psi,
            const structure::QuantumSystem<RANK>& sys,
            const ParsEvolution& p)
{
  evolve<V>(psi,cpputils::sharedPointerize(sys),p);
}


template<int RANK>
inline
void evolve(quantumdata::StateVector<RANK>& psi,
            const structure::QuantumSystem<RANK>& sys,
            const ParsEvolution& p)
{
  evolve<tmptools::V_Empty>(psi,cpputils::sharedPointerize(sys),p);
}


template<int V0, int... V_REST, typename SV, typename SYS>
inline
void evolve(SV& psi,
            const SYS& sys,
            const ParsEvolution& p)
{
  evolve<tmptools::Vector<V0,V_REST...> >(psi,sys,p);
}


#endif // QUANTUMTRAJECTORY_EVOLUTION__H_INCLUDED
