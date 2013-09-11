// -*- C++ -*-
/// \briefFile{Forward-declares the general form of structure::Exact}
#ifndef   STRUCTURE_EXACTFWD_H_INCLUDED
#define   STRUCTURE_EXACTFWD_H_INCLUDED

namespace structure {


class ExactCommon;

/// The interface every system that needs transformation between two quantum mechanical pictures must present towards the trajectory drivers
/**
 * Experience shows that even when a system uses interaction picture (which is automatically the case if any of its subsystems does) – that is, part of its dynamics is solved exactly – 
 * it may still want to calculate the jump operators and quantum averages in the normal picture. (cf. Cases 1 & 3 \link TimeDependenceLevel above\endlink.)
 * This is useful e.g. to reuse the code written for the non-interaction-picture case.
 * 
 * In this case, the framework has to be provided with some means to transform between the two pictures.
 * This is fulfilled by this class, from which classes describing such systems have to inherit.
 * 
 * E.g. if quantumtrajectory::MCWF_Trajectory sees that the simulated system inherits from Exact, then it will make the coherent part of the evolution in interaction picture,
 * whereupon it transforms back to normal picture, so that all the rest (jump rates, eventual jumps, quantum averages) can be calculated in this latter picture.
 * This makes that the two pictures coincide before each timestep. (Cf. also the stages described @ quantumtrajectory::MCWF_trajectory.)
 * 
 * The design is very similar to that of Hamiltonian, the general template never being defined.
 * 
 * \tparamRANK
 * \tparam IS_TWO_TIME default `true`, the most general case
 * 
 * \see Exact<RANK,true>, Exact<RANK,false> in Exact.h
 * 
 */
template<int, bool IS_TWO_TIME=true>
class Exact;


} // structure

#endif // STRUCTURE_EXACTFWD_H_INCLUDED
