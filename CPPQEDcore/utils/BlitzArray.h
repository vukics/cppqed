// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines template aliases for real and complex arrays + implements the traits functions declared in ArrayTraits.h for `blitz::Array`}
/** \todo Much of the traits in this file are superfluous as many components can be formulated generally for all `blitz::Array` types */
#ifndef   CPPQEDCORE_UTILS_BLITZARRAY_H_INCLUDED
#define   CPPQEDCORE_UTILS_BLITZARRAY_H_INCLUDED

#include "ArrayTraits.h"
#include "ComplexExtensions.h"
#include "SliceIterator.h"

#include <blitz/array.h>


/// An array of doubles of arbitrary arity
template <int RANK> using DArray=blitz::Array<double,RANK>;

/// A complex array of arbitrary arity
template <int RANK> using CArray=blitz::Array<dcomp ,RANK>;


namespace cpputils {

/** \cond SPECIALIZATION */

namespace sliceiterator {


template<typename T, int RANK, typename ... SubscriptPack>
auto subscript(const blitz::Array<T,RANK>& array, SubscriptPack&&... subscriptPack)
{
  static_assert( sizeof...(SubscriptPack)==RANK , "Incorrect number of subscripts for blitz::Array." );
  return array(std::forward<SubscriptPack>(subscriptPack)...);
}


/**
 * \todo SliceIterator could be defined as taking an ARRAY *type* as template param instead of a *template* (a traits metafunction calculating ResArray in this case)
 * Then, the begin, end, fullRange functions could simply take also an ARRAY template param, that would greatly ease template param deduction,
 * making also the following specializations unnecessary.
 */
template<typename V, int RANK>
auto begin(const CArray<RANK>& array) {return begin<V,CArray>(array);}

template<typename V, int RANK>
auto end  (const CArray<RANK>& array) {return end<V,CArray>(array);}

template<typename V, int RANK>
auto fullRange(const CArray<RANK>& array) {return fullRange<V,CArray>(array);}


} // sliceiterator



template<typename Numtype, int RANK>
struct Rank<blitz::Array<Numtype,RANK> > : boost::mpl::int_<RANK> {};

template<int RANK>
struct TypeID<DArray<RANK> >
{
  static const std::string value;
};

template<int RANK>
const std::string TypeID<DArray<RANK> >::value="DArray";


template<int RANK>
struct TypeID<CArray<RANK> >
{
  static const std::string value;
};

template<int RANK>
const std::string TypeID<CArray<RANK> >::value="CArray";

/** \endcond */


/// \name `blitz::Array` memory traits for `blitz::Array<double,n>`
//@{

template<int n>
inline bool isStorageContiguous(const DArray<n>& a) {return a.isStorageContiguous();}


template<int n>
inline size_t size(const DArray<n>& a) {return a.size();}


template<int n>
inline std::vector<size_t> dimensions(const DArray<n>& a) {return std::vector<size_t>(a.extent().begin(),a.extent().end());}


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


template<int n>
inline std::vector<size_t> dimensions(const CArray<n>& a) {return std::vector<size_t>(a.extent().begin(),a.extent().end());}


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


#endif // CPPQEDCORE_UTILS_BLITZARRAY_H_INCLUDED
