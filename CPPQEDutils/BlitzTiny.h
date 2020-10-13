// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines template aliases for `blitz::TinyVector`s used for characterising the size of multi-arrays and indexing them}
#ifndef   CPPQEDCORE_UTILS_BLITZTINY_H_INCLUDED
#define   CPPQEDCORE_UTILS_BLITZTINY_H_INCLUDED

#include <blitz/array.h>

/// A tiny vector describing extensions of objects of arbitrary arity
template <int RANK> using ExtTiny=blitz::TinyVector<   size_t,RANK>;

/// A tiny vector used for indexing of objects of arbitrary arity
template <int RANK> using IdxTiny=blitz::TinyVector<ptrdiff_t,RANK>;

#endif // CPPQEDCORE_UTILS_BLITZTINY_H_INCLUDED
