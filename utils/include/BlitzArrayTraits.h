// -*- C++ -*-
#ifndef UTILS_INCLUDE_BLITZARRAYTRAITS_H_INCLUDED
#define UTILS_INCLUDE_BLITZARRAYTRAITS_H_INCLUDED

#include "BlitzArrayTraitsFwd.h"


namespace blitzplusplus {

#define RETURN_type blitz::Array<T,RANK>

namespace details {

template<typename T, int RANK>
class MyArray : public RETURN_type
{
public:
  static const blitz::GeneralArrayStorage<RANK>& stealStorage(const RETURN_type& array) {return array.storage_;}
};

} // details

template<typename T, int RANK>
const RETURN_type
clone(T* restrict data, const RETURN_type& array)
{
  RETURN_type(data,array.shape(),array.stride(),blitz::neverDeleteData/*,details::MyArray<T,RANK>::stealStorage(array)*/);
}

#undef  RETURN_type

} // blitzplusplus


/////////////////////////
//
// blitz::Array<double,n>
//
/////////////////////////

namespace cpputils {

template<int n>
struct ArrayMemoryTraits<TTD_DARRAY(n)> {

  typedef TTD_DARRAY(n) DA;
  
  static bool isStorageContiguous(const DA& a) {return a.isStorageContiguous();}

  static size_t size(const DA& a) {return a.size();}

  static const double* data(const DA& a) {return a.size() ? a.data() : 0;}
  static       double* data(      DA& a) {return const_cast<double*>(data(static_cast<const DA&>(a)));}

  static       DA create(      double* y, const DA& a) {return DA(y,a.shape(),blitz::neverDeleteData);} // blitzplusplus::clone(y,a);}
  static const DA create(const double* y, const DA& a) {return create(const_cast<double*>(y),a);}

  static DA create(const DA& a) {return DA(a.shape());}

};


////////////////////////
//
// blitz::Array<dcomp,n>
//
////////////////////////


template<int n>
struct ArrayMemoryTraits<TTD_CARRAY(n)> {

  typedef TTD_CARRAY(n) CA;
  
  static bool isStorageContiguous(const CA& a) {return a.isStorageContiguous();}

  static size_t size(const CA& a) {return a.size()<<1;} // The size of the underlying double* storage!!!

  // static const double* data(const CA& a) {return a.size() ? reinterpret_cast<const double*>(a.data()) : 0;}
  static const double* data(const CA& a) {return a.size() ? real(a).data() : 0;}
  static       double* data(      CA& a) {return const_cast<double*>(data(static_cast<const CA&>(a)));}

  static       CA create(      double* y, const CA& a) {return CA(reinterpret_cast<dcomp*>(y),a.shape(),blitz::neverDeleteData);}
  // blitzplusplus::clone(reinterpret_cast<dcomp*>(y),a);}
  static const CA create(const double* y, const CA& a) {return create(const_cast<double*>(y),a);}

  static CA create(const CA& a) {return CA(a.shape());}

};




template<>
struct ArrayTraversalTraits<TTD_DARRAY(1)> : ArrayMemoryTraits<TTD_DARRAY(1)> {

  typedef TTD_DARRAY(1) DA;

  static const double& ss(const DA& a, size_t i) {return a(i);}
  static       double& ss(      DA& a, size_t i) {return const_cast<double&>(ss(static_cast<const DA&>(a),i));}

  static size_t ssLimit(const DA& a) {return a.size();}

};


template<int n>
struct ArrayTraversalTraits<TTD_CARRAY(n)> : ArrayMemoryTraits<TTD_CARRAY(n)> {

  typedef TTD_CARRAY(n) CA;

  static const dcomp& ss(const CA& a, size_t i) {return *(a.begin()+i);}
  // This funny-looking solution is aimed at ensuring that the array is traversed contiguously independently of the rank. Bounds checking? --- CheckedIterator.
  static       dcomp& ss(      CA& a, size_t i) {return const_cast<dcomp&>(ss(static_cast<const CA&>(a),i));}

  static size_t ssLimit(const CA& a) {return a.size();}

};


template<> 
struct ArrayTraversalTraits<TTD_CARRAY(1)> : ArrayMemoryTraits<TTD_CARRAY(1)> {

  typedef TTD_CARRAY(1) CA;

  static const dcomp& ss(const CA& a, size_t i) {return a(i);}
  static       dcomp& ss(      CA& a, size_t i) {return const_cast<dcomp&>(ss(static_cast<const CA&>(a),i));}

  static size_t ssLimit(const CA& a) {return a.size();}

  static size_t stride(const CA& a) {return a.stride(0);}

};

} // cpputils

#endif // UTILS_INCLUDE_BLITZARRAYTRAITS_H_INCLUDED
