// -*- C++ -*-
/// \briefFile{Defines enumeration structure::TimeDependence & forward-declares the general form of structure::Hamiltonian}
#ifndef   STRUCTURE_HAMILTONIANFWD_H_INCLUDED
#define   STRUCTURE_HAMILTONIANFWD_H_INCLUDED

namespace structure {

/// Enumeration of different possibilities for time dependence of Hamiltonians
/** With \f$t_0\f$ being the time instant where the \link Exact two pictures\endlink coincide: */
enum TimeDependence {
  TWO_TIME, ///< **Case 1** \f$H(t,t_0)\f$ – Time-dependent problem + \link Exact exact part\endlink (\f$U(t,t_0)\f$)
  ONE_TIME, ///< **Case 2** \f$H(t)\f$ – Time-dependent problem, no exact part **OR Case 3** \f$H(t-t_0)\f$ – Time-independent problem + \link Exact exact part\endlink (\f$U(t-t_0)\f$)
  NO_TIME   ///< **Case 4** \f$H(0)\f$ – Time-independent problem, no exact part
};


/// The interface every system having (possibly non-Hermitian) Hamiltonian time-evolution must present towards the trajectory drivers
/**
 * The most general template is never defined, only its partial specializations in the second parameter.
 * 
 * \tparamRANK
 * \tparam TD Degree of \link TimeDependence time dependence\endlink. The most general (#TWO_TIME) is taken as default.
 * 
 * \see Hamiltonian<RANK,TWO_TIME> through Hamiltonian<RANK,NO_TIME> in Hamiltonian.h
 * 
 */
template<int RANK,TimeDependence TD=TWO_TIME>
class Hamiltonian;


} // structure

#endif // STRUCTURE_HAMILTONIANFWD_H_INCLUDED
