/// \briefFile{Declarations of traits functions for adapting array types to generic functions}
// -*- C++ -*-
#ifndef   UTILS_ARRAYTRAITS_H_INCLUDED
#define   UTILS_ARRAYTRAITS_H_INCLUDED

#include <vector>

namespace cpputils {



/// \name Array memory traits
//@{


template<typename A>
bool isStorageContiguous(const A& a);


template<typename A>
size_t size(const A& a);

template<typename A>
std::vector<size_t> dimensions(const A& a);

template<typename A>
size_t rank();


template<typename A>
inline size_t rank(const A& a) {return rank<A>();};



template<typename A>
const double* data(const A& a);

template<typename A>
      double* data(      A& a);



template<typename A>
      A create(      double* y, const A& a); ///< Clone (create a non-owning array of data `y` of the same memory layout as `a`)

template<typename A>
const A create(const double* y, const A& a); ///< Const clone (create a const non-owning array of data `y` of the same memory layout as `a`)


template<typename A>
A create(const A& a); ///< Empty clone (create a newly allocated owning empty array of the same memory layout as `a`)
//@}



/// \name Array traversal traits
//@{

/// subscription of `a` (which might be a multi-array) with a *single* integer
/**
 * \note `A::element_type` is a hypothetic member type, which in fact never plays a role, since we are relying on template-parameter inference & overload resolution,
 * wherein the return type does not play any role
 */
template<typename A>
const typename A::element_type& subscript(const A& a, size_t i);


/// non-const subscription
template<typename A>
      typename A::element_type& subscript(      A& a, size_t i);


template<typename A>
size_t subscriptLimit(const A& a);
  
//@}
  
} // cpputils


#endif // UTILS_ARRAYTRAITS_H_INCLUDED
