// -*- C++ -*-
#ifndef   UTILS_INCLUDE_ARRAYTRAITS_H_INCLUDED
#define   UTILS_INCLUDE_ARRAYTRAITS_H_INCLUDED

namespace cpputils {


/////////////////////////
//
// „Memory traits”
//
/////////////////////////


template<typename A>
bool isStorageContiguous(const A& a);


template<typename A>
size_t size(const A& a);


template<typename A>
const double* data(const A& a);

template<typename A>
      double* data(      A& a);


template<typename A>
      A create(      double* y, const A& a);

template<typename A>
const A create(const double* y, const A& a);


template<typename A>
A create(const A& a);


/////////////////////////
//
// „Traversal traits”
//
/////////////////////////

template<typename A>
const typename A::element_type& subscript(const A& a, size_t i);
// A::element_type is a hypothetic member type, which in fact never plays a role, since we are relying on template parameter inference & overload resolution

template<typename A>
      typename A::element_type& subscript(      A& a, size_t i);


template<typename A>
size_t subscriptLimit(const A& a);
  
  
} // cpputils


#endif // UTILS_INCLUDE_ARRAYTRAITS_H_INCLUDED
