// -*- C++ -*-
#ifndef   _BLITZ_ARRAY_TRAITS_FWD_H
#define   _BLITZ_ARRAY_TRAITS_FWD_H

#include "ArrayTraitsFwd.h"
#include "BlitzArray.h"


namespace cpputils {

template<int n>
struct ArrayMemoryTraits<TTD_DARRAY(n)>;

template<int n>
struct ArrayMemoryTraits<TTD_CARRAY(n)>;

template<>
struct ArrayTraversalTraits<TTD_DARRAY(1)>;

template<int n>
struct ArrayTraversalTraits<TTD_CARRAY(n)>;

template<> 
struct ArrayTraversalTraits<TTD_CARRAY(1)>;

} // cpputils

#endif // _BLITZ_ARRAY_TRAITS_FWD_H
