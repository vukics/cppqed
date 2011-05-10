// -*- C++ -*-
#ifndef   STRUCTURE_HAMILTONIAN_FWD_INCLUDED
#define   STRUCTURE_HAMILTONIAN_FWD_INCLUDED

namespace structure {


enum TimeDependence {TWO_TIME, ONE_TIME, NO_TIME};


template<int,TimeDependence=TWO_TIME>
class Hamiltonian;


} // structure

#endif // STRUCTURE_HAMILTONIAN_FWD_INCLUDED
