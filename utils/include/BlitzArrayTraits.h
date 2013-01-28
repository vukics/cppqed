// -*- C++ -*-
#ifndef UTILS_INCLUDE_BLITZARRAYTRAITS_H_INCLUDED
#define UTILS_INCLUDE_BLITZARRAYTRAITS_H_INCLUDED

#include "BlitzArray.h"


namespace cpputils {

/////////////////////////
//
// „Memory traits”
//
/////////////////////////


// blitz::Array<double,n>


template<int n>
inline bool isStorageContiguous(const TTD_DARRAY(n)& a) {return a.isStorageContiguous();}


template<int n>
inline size_t size(const TTD_DARRAY(n)& a) {return a.size();}


template<int n>
inline const double* data(const TTD_DARRAY(n)& a) {return a.size() ? a.data() : 0;}

template<int n>
inline       double* data(      TTD_DARRAY(n)& a) {return const_cast<double*>(data(static_cast<const TTD_DARRAY(n)&>(a)));}


template<int n>
inline       TTD_DARRAY(n) create(      double* y, const TTD_DARRAY(n)& a) {return TTD_DARRAY(n)(y,a.shape(),blitz::neverDeleteData);}

template<int n>
inline const TTD_DARRAY(n) create(const double* y, const TTD_DARRAY(n)& a) {return create(const_cast<double*>(y),a);}


template<int n>
inline TTD_DARRAY(n) create(const TTD_DARRAY(n)& a) {return TTD_DARRAY(n)(a.shape());}


// blitz::Array<dcomp,n>


template<int n>
inline bool isStorageContiguous(const TTD_CARRAY(n)& a) {return a.isStorageContiguous();}


template<int n>
inline size_t size(const TTD_CARRAY(n)& a) {return a.size()<<1;} // The size of the underlying double* storage!!!


template<int n>
inline const double* data(const TTD_CARRAY(n)& a) {return a.size() ? real(a).data() : 0;}

template<int n>
inline       double* data(      TTD_CARRAY(n)& a) {return const_cast<double*>(data(static_cast<const TTD_CARRAY(n)&>(a)));}


template<int n>
inline       TTD_CARRAY(n) create(      double* y, const TTD_CARRAY(n)& a) {return TTD_CARRAY(n)(reinterpret_cast<dcomp*>(y),a.shape(),blitz::neverDeleteData);}

template<int n>
inline const TTD_CARRAY(n) create(const double* y, const TTD_CARRAY(n)& a) {return create(const_cast<double*>(y),a);}


template<int n>
inline TTD_CARRAY(n) create(const TTD_CARRAY(n)& a) {return TTD_CARRAY(n)(a.shape());}



/////////////////////////
//
// „Traversal traits”
//
/////////////////////////


inline const double& subscript(const TTD_DARRAY(1)& a, size_t i) {return a(i);}
inline       double& subscript(      TTD_DARRAY(1)& a, size_t i) {return const_cast<double&>(subscript(static_cast<const TTD_DARRAY(1)&>(a),i));}

inline size_t subscriptLimit(const TTD_DARRAY(1)& a) {return a.size();}


template<int n>
inline const dcomp& subscript(const TTD_CARRAY(n)& a, size_t i) {return *(a.begin()+i);}
  // This funny-looking solution is aimed at ensuring that the array is traversed contiguously independently of the rank. Bounds checking? --- CheckedIterator.

template<int n>
inline       dcomp& subscript(      TTD_CARRAY(n)& a, size_t i) {return const_cast<dcomp&>(subscript(static_cast<const TTD_CARRAY(n)&>(a),i));}

template<int n>
inline size_t subscriptLimit(const TTD_CARRAY(n)& a) {return a.size();}


inline const dcomp& subscript(const TTD_CARRAY(1)& a, size_t i) {return a(i);}
inline       dcomp& subscript(      TTD_CARRAY(1)& a, size_t i) {return const_cast<dcomp&>(subscript(static_cast<const TTD_CARRAY(1)&>(a),i));}

inline size_t subscriptLimit(const TTD_CARRAY(1)& a) {return a.size();}

inline size_t stride(const TTD_CARRAY(1)& a) {return a.stride(0);}


} // cpputils

#endif // UTILS_INCLUDE_BLITZARRAYTRAITS_H_INCLUDED
