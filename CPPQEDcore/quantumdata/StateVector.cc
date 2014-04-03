// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

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
