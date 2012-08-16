// -*- C++ -*-
#ifndef   UTILS_INCLUDE_BLITZARRAYTRAITSFWD_H_INCLUDED
#define   UTILS_INCLUDE_BLITZARRAYTRAITSFWD_H_INCLUDED

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

#endif // UTILS_INCLUDE_BLITZARRAYTRAITSFWD_H_INCLUDED
