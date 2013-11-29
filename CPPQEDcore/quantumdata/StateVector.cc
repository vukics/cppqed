
// explicitely instantiate debug function if in debug mode

#ifndef NDEBUG
#include "StateVector.tcc"

namespace quantumdata {

template class StateVector<1>;
template class StateVector<2>;
template class StateVector<3>;
template class StateVector<4>;

}
#endif // NDEBUG
