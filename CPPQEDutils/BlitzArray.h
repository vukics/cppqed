// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines template aliases for real and complex arrays + implements the traits functions declared in ArrayTraits.h for `blitz::Array`}
#ifndef   CPPQEDCORE_UTILS_BLITZARRAY_H_INCLUDED
#define   CPPQEDCORE_UTILS_BLITZARRAY_H_INCLUDED

#include "ArrayTraits.h"
#include "ComplexExtensions.h"
#include "SliceIterator.h"

#include <boost/numeric/odeint.hpp>

#include <blitz/array.h>



namespace boost { namespace numeric { namespace odeint {
  
template<typename Numtype, int RANK>
struct is_resizeable<blitz::Array<Numtype,RANK>> : boost::true_type {};

template<typename Numtype, int RANK>
struct same_size_impl<blitz::Array<Numtype,RANK>, blitz::Array<Numtype,RANK>>
{ // define how to check size
  static bool same_size(const blitz::Array<Numtype,RANK> &v1, const blitz::Array<Numtype,RANK> &v2) {return all( v1.shape() == v2.shape() );}
};

template<typename Numtype, int RANK>
struct resize_impl<blitz::Array<Numtype,RANK>, blitz::Array<Numtype,RANK>>
{ // define how to resize
  static void resize(blitz::Array<Numtype,RANK> &v1, const blitz::Array<Numtype,RANK> &v2) {v1.resize( v2.shape() );}
};

template<typename Numtype, int RANK>
struct vector_space_norm_inf<blitz::Array<Numtype,RANK>>
{
  typedef double result_type;
  double operator()(const blitz::Array<Numtype,RANK>& v ) const
  {
    return max( abs(v) );
  }
};

template<typename Numtype, int RANK>
struct norm_result_type<blitz::Array<Numtype,RANK>> : mpl::identity<double> {};

} } } // boost::numeric::odeint


/// An array of doubles of arbitrary arity
template <int RANK> using DArray=blitz::Array<double,RANK>;

/// A complex array of arbitrary arity
template <int RANK> using CArray=blitz::Array<dcomp ,RANK>;


namespace cppqedutils {

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


template<typename T, int RANK>
struct ElementType<blitz::Array<T,RANK>> : boost::mpl::identity<T> {};

template<typename Numtype, int RANK>
constexpr int Rank_v<blitz::Array<Numtype,RANK>> = RANK;

template<int RANK>
constexpr auto TypeID_v<DArray<RANK>> ="DArray";

template<int RANK>
constexpr auto TypeID_v<CArray<RANK>> ="CArray";

/** \endcond */


template<typename T, int n>
struct IsStorageContiguous<blitz::Array<T,n>>
{
  static auto _(const blitz::Array<T,n>& a) {return a.isStorageContiguous();}
};


template<typename T, int n>
struct Extents<blitz::Array<T,n>>
{
  static auto _(const blitz::Array<T,n>& a) {Extents_t<n> res; std::copy(a.extent().begin(),a.extent().end(),res.begin()); return res;}
};


/// \name `blitz::Array` memory traits for `blitz::Array<double,n>`
//@{

template<int n>
struct Data_c<DArray<n>>
{
  static const double* _(const DArray<n>& a) {return a.size() ? a.data() : nullptr;}
};


template<int n>
struct Create<DArray<n>>
{
  static DArray<n> _(double* y, Extents_t<n> e) {ExtTiny<n> et; std::copy(e.begin(),e.end(),et.begin()); return DArray<n>(y,et,blitz::neverDeleteData);}
};


//@}


/// \name `blitz::Array` memory traits for `blitz::Array<dcomp,n>`
//@{

template<int n>
struct Size<CArray<n>>
{
  static auto _(const CArray<n>& a) {return a.size()<<1;} // The size of the underlying double* storage!!!
};


template<int n>
struct Data_c<CArray<n>>
{
  static const double* _(const CArray<n>& a) {return a.size() ? real(a).data() : nullptr;}
};


template<int n>
struct Create<CArray<n>>
{
  static CArray<n> _(double* y, Extents_t<n> e)
  {
    ExtTiny<n> et; std::copy(e.begin(),e.end(),et.begin());
    return CArray<n>(reinterpret_cast<dcomp*>(y),et,blitz::neverDeleteData);    
  }
};


//@}

template<typename T, int n>
struct CreateFromExtents<blitz::Array<T,n>>
{
  static blitz::Array<T,n> _(Extents_t<n> e) {ExtTiny<n> et; std::copy(e.begin(),e.end(),et.begin()); return blitz::Array<T,n>(et);}
};


template<typename T, int n>
struct Copy<blitz::Array<T,n>>
{
  /// Copy of a with its own memory
  static auto _(const blitz::Array<T,n>& a) {return a.copy();}
};


/// \name `blitz::Array` traversal traits for unary double and complex arrays
//@{

template<typename T>
struct Subscript_c<blitz::Array<T,1>>
{
  static const auto& _(const blitz::Array<T,1>& a, size_t i) {return a(i);}
};


template<typename T>
struct Stride<blitz::Array<T,1>>
{
  static auto _(const blitz::Array<T,1>& a) {return a.stride(0);}
};

//@}



/// \name `blitz::Array` traversal traits for double and complex arrays of arbitrary arity
//@{
/**
 * \note this is broken since `blitz::Array` iterators are not random-access
 * \todo fix this
 */
template<int n>
inline const dcomp& subscript(const CArray<n>& a, size_t i) {return *(a.begin()+i);}

template<int n>
inline size_t subscriptLimit(const CArray<n>& a) {return a.size();}
//@}

} // cppqedutils


#endif // CPPQEDCORE_UTILS_BLITZARRAY_H_INCLUDED
