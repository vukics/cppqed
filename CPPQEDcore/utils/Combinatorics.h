// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDCORE_UTILS_COMBINATORICS_H_INCLUDED
#define   CPPQEDCORE_UTILS_COMBINATORICS_H_INCLUDED

#include "CombinatoricsFwd.h"

#include "Exception.h"

#include <blitz/array.h>

#include <boost/range/algorithm/equal.hpp>

namespace cpputils {


class CWR_Dir // Combinations with repetitions
{
public:
  typedef blitz::Array<size_t,2> Impl;
  typedef blitz::Array<size_t,1> Configuration; 

  CWR_Dir(size_t n, size_t k);
  // n is the number of objects from which you can choose and k is the number to be chosen

  // E.g. for bosonic states on a lattice: n is the number of lattice sites, k is the number of bosonic particles.

  const Configuration operator[](size_t  ) const;
  const Configuration operator[](int    i) const {return operator[](size_t(i));}
  
  enum SubscriptingProblem {SIZE_MISMATCH, NOT_FOUND};

  template<typename C, SubscriptingProblem P>
  struct SubscriptingException : public Exception
  {
    SubscriptingException(const C& c) : conf(c) {}
    
    const C conf;
  };

  template<typename C>
  size_t operator[](const C& c) const
  // gives the ordinal number of a certain configuration, which is given by a Boost.Range-complying range C
  {
    if (c.size()!=impl_.extent(1)) throw SubscriptingException<C,SIZE_MISMATCH>(c);
    for (size_t i=0; i<impl_.extent(0); i++)
      if (boost::equal(c,impl_(i,blitz::Range::all()))) return i;
    throw SubscriptingException<C,NOT_FOUND>(c);
  }

  const Impl& operator()() const {return impl_;}

private:
  static const Impl recurse(Impl, size_t);

  const Impl impl_;

};




} // cpputils


#endif // CPPQEDCORE_UTILS_COMBINATORICS_H_INCLUDED
