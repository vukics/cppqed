// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Definition of cpputils::SliceIterator together with its helpers}
#if !BOOST_PP_IS_ITERATING

#ifndef CPPQEDCORE_UTILS_SLICEITERATOR_TCC_INCLUDED
#define CPPQEDCORE_UTILS_SLICEITERATOR_TCC_INCLUDED

#include "SliceIterator.h"

#include <boost/preprocessor/repetition/enum.hpp>



namespace cpputils { namespace sliceiterator {

//////////////////////////
//
// Indexer Implementations
//
//////////////////////////

namespace details {

// This class serves as private base for Indexer
template<int RANK, typename V>
class IndexerBase
{
public:
  typedef typename boost::mpl::fold<tmptools::Ordinals<RANK>,
                                    boost::fusion::vector<>,
                                    boost::mpl::push_back<boost::mpl::_1,
                                                          boost::mpl::if_<tmptools::numerical_contains<V,boost::mpl::_2>,blitz::Range,int>
                                                          >
                                    >::type Idx;

private:
  typedef typename VecIdxTiny<RANK,V>::const_iterator CI;

public:
  static const Idx&
  fillIdxValues(const VecIdxTiny<RANK,V>& idx)
  {
    boost::fusion::fold(cache_,idx.begin(),
                        [](CI iter, auto& t) -> CI {
                          using T=std::decay_t<decltype(t)>;
                          if constexpr (std::is_same<T,blitz::Range>::value) {t=blitz::Range::all(); return iter;}
                          else if constexpr (std::is_same<T,int>::value) {t=*iter++; return iter;}
                        });
    return cache_;
  }

protected:
  static Idx cache_; 
  // In this way the creation of a default-constructed Helper & Idx objects in fillIdxValues can be avoided

};

// Definition of members cache_ (the above is only a declaration!):

template<int RANK, typename V>
typename IndexerBase<RANK,V>::Idx IndexerBase<RANK,V>::cache_;


////////////////////////////
//
// Transpose Implementations
//
////////////////////////////


namespace namehider {

using namespace boost::mpl;
namespace mpl=boost::mpl;
using namespace tmptools;

template<int RANK, typename V>
struct Algorithm 
  : fold<Ordinals<RANK>,
         pair<vector_c<int>,typename boost::mpl::begin<V>::type>,
         pair<push_back<mpl::first<mpl::_1>,
                        if_<numerical_contains<V,mpl::_2>,
                            deref<second<mpl::_1> >,
                            mpl::_2
                            >
                        >,
              if_<numerical_contains<V,mpl::_2>,
                  mpl::next<second<mpl::_1> >,
                  second<mpl::_1>
                  >
              >
         >
{};

} // namehider


template<int RANK, typename V>
struct TransposerMeta : boost::mpl::first<typename namehider::Algorithm<RANK,V>::type>
{};

} // details


} } // cpputils::sliceiterator



//////////////////////////////////////////
//
// Brute Force Implementations for Indexer
//
//////////////////////////////////////////


#define BOOST_PP_ITERATION_LIMITS (1,BLITZ_ARRAY_LARGEST_RANK)
#define BOOST_PP_FILENAME_1 "SliceIterator.tcc"

#include BOOST_PP_ITERATE()

#undef BOOST_PP_FILENAME_1
#undef BOOST_PP_ITERATION_LIMITS


#endif // CPPQEDCORE_UTILS_SLICEITERATOR_TCC_INCLUDED


#else  // BOOST_PP_IS_ITERATING


#define rank BOOST_PP_ITERATION()

#define INDEXER_print1(z,m,data) boost::mpl::at_c<TM,m>::type::value
#define INDEXER_print2(z,m,data) boost::fusion::at_c<m>(Base::cache_)

namespace cpputils { namespace sliceiterator {

template<template <int> class ARRAY, typename V> struct Transposer<ARRAY,rank,V>
{
  typedef typename details::TransposerMeta<rank,V>::type TM;

  static auto& _(ARRAY<rank>& array)
  {
    array.transposeSelf(BOOST_PP_ENUM(rank,INDEXER_print1,~) );
    return array;
  }

};


template<template <int> class ARRAY, typename V> struct Indexer<ARRAY,rank,V> : private details::IndexerBase<rank,V>
{
  typedef details::IndexerBase<rank,V> Base;

  static auto& _(const ARRAY<rank>& array, ResArray<ARRAY,V>& resArray, const VecIdxTiny<rank,V>& idx)
  {
    Base::fillIdxValues(idx);
    resArray.reference(subscript(array,BOOST_PP_ENUM(rank,INDEXER_print2,~) ) );
    return resArray;
  }

};

} } // cpputils::sliceiterator

#undef INDEXER_print1
#undef INDEXER_print2

#undef rank


#endif // BOOST_PP_IS_ITERATING
