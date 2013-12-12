/// \briefFile{Implementations of the traits functions declared in ArrayTraits.h for `blitz::Array`}
// -*- C++ -*-
#ifndef UTILS_BLITZARRAYTRAITS_H_INCLUDED
#define UTILS_BLITZARRAYTRAITS_H_INCLUDED

#include "BlitzArray.h"
#include <boost/iterator/iterator_concepts.hpp>


namespace cpputils {


/// \name `blitz::Array` memory traits for `blitz::Array<double,n>`
//@{


template<int n>
inline bool isStorageContiguous(const DArray<n>& a) {return a.isStorageContiguous();}


template<int n>
inline size_t size(const DArray<n>& a) {return a.size();}


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
template<int n>
inline size_t rank(const DArray<n>& a) {return n;} // The parameter a is actually only needed for template argument deduction.
#pragma GCC diagnostic pop


template<int n>
inline const double* data(const DArray<n>& a) {return a.size() ? a.data() : 0;}

template<int n>
inline       double* data(      DArray<n>& a) {return const_cast<double*>(data(static_cast<const DArray<n>&>(a)));}


template<int n>
inline       DArray<n> create(      double* y, const DArray<n>& a) {return DArray<n>(y,a.shape(),blitz::neverDeleteData);}

template<int n>
inline const DArray<n> create(const double* y, const DArray<n>& a) {return create(const_cast<double*>(y),a);}


template<int n>
inline DArray<n> create(const DArray<n>& a) {return DArray<n>(a.shape());}
//@}


/// \name `blitz::Array` memory traits for `blitz::Array<dcomp,n>`
//@{

template<int n>
inline bool isStorageContiguous(const CArray<n>& a) {return a.isStorageContiguous();}


template<int n>
inline size_t size(const CArray<n>& a) {return a.size()<<1;} // The size of the underlying double* storage!!!


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
template<int n>
inline size_t rank(const CArray<n>& a) {return n;} // The parameter a is actually only needed for template argument deduction.
#pragma GCC diagnostic pop


template<int n>
inline const double* data(const CArray<n>& a) {return a.size() ? real(a).data() : 0;}

template<int n>
inline       double* data(      CArray<n>& a) {return const_cast<double*>(data(static_cast<const CArray<n>&>(a)));}


template<int n>
inline       CArray<n> create(      double* y, const CArray<n>& a) {return CArray<n>(reinterpret_cast<dcomp*>(y),a.shape(),blitz::neverDeleteData);}

template<int n>
inline const CArray<n> create(const double* y, const CArray<n>& a) {return create(const_cast<double*>(y),a);}


template<int n>
inline CArray<n> create(const CArray<n>& a) {return CArray<n>(a.shape());}
//@}



/// \name `blitz::Array` traversal traits for unary double and complex arrays
//@{

inline const double& subscript(const DArray<1>& a, size_t i) {return a(i);}
inline       double& subscript(      DArray<1>& a, size_t i) {return const_cast<double&>(subscript(static_cast<const DArray<1>&>(a),i));}

inline size_t subscriptLimit(const DArray<1>& a) {return a.size();}


inline const dcomp& subscript(const CArray<1>& a, size_t i) {return a(i);}
inline       dcomp& subscript(      CArray<1>& a, size_t i) {return const_cast<dcomp&>(subscript(static_cast<const CArray<1>&>(a),i));}

inline size_t subscriptLimit(const CArray<1>& a) {return a.size();}

inline size_t stride(const CArray<1>& a) {return a.stride(0);}
//@}



/// \name `blitz::Array` traversal traits for unary double and complex arrays
//@{
/**
 * \note this is broken since `blitz::Array` iterators are not random-access
 * \todo fix this
 */
template<int n>
inline const dcomp& subscript(const CArray<n>& a, size_t i) {return *(a.begin()+i);}

template<int n>
inline       dcomp& subscript(      CArray<n>& a, size_t i) {return const_cast<dcomp&>(subscript(static_cast<const CArray<n>&>(a),i));}

template<int n>
inline size_t subscriptLimit(const CArray<n>& a) {return a.size();}
//@}

} // cpputils

#endif // UTILS_BLITZARRAYTRAITS_H_INCLUDED
