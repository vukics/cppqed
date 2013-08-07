// -*- C++ -*-
/// \briefFile{Forward-declares the general structure::Liouvillean}
#ifndef   STRUCTURE_LIOUVILLEANFWD_H_INCLUDED
#define   STRUCTURE_LIOUVILLEANFWD_H_INCLUDED


namespace structure {


/// The interface every system having Liouvillean time-evolution must present towards the trajectory drivers
/**
 * The time-evolution must be Markovian where the Lindblad form of the Master equation is the most general one possible:
 * \f[\dot\rho=2\real{\frac\HnH{i\hbar}\rho}+\sum_mJ_m\rho J_m^\dagger\f]
 * The class represents the set of \f$J_m\f$ operators (Lindblads or quantum jump operators) of arbitrary number,
 * either for Master equation or Monte Carlo wave-function evolution (in the latter case it calculates the jump rates as well).
 * 
 * \tparamRANK
 * \tparam IS_TIME_DEPENDENT describes whether the \f$J_m\f$s are time-dependent. Default `true`, the most general case.
 * 
 * \note No matter how complex the quantum system, the framework always assigns a unique ordinal to each jump corresponding to every subsystem
 * 
 * Similarly to Hamiltonian & Exact, the general template is never defined, but a hierarchy of partial specializations in the second template argument.
 * 
 * \see Liouvillean<RANK,true>, Liouvillean<RANK,false> in Liouvillean.h
 * 
 * \note It is always possible to forgo the explicit calculation of certain jump rates because the probability can be calculated also on the basis of  the actWithJ() by the MCWF stepper. The fact that such a fallback is desired can be signalled by setting a negative value for the probability of the given jump (“special jump”).
 * 
 */
template<int RANK, bool IS_TIME_DEPENDENT=true>
class Liouvillean;


} // structure


#endif // STRUCTURE_LIOUVILLEANFWD_H_INCLUDED
