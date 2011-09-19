// -*- C++ -*-
#ifndef   TMP_HELPERS_INCLUDED
#define   TMP_HELPERS_INCLUDED

#include <boost/mpl/max_element.hpp>


namespace composite {

using namespace mpl;

template<typename VEC, typename Pred=less<mpl::_1,mpl::_2> >
struct MaxMF : deref<typename max_element<VEC,Pred>::type>::type {};

template<typename VEC1, typename VEC2>
struct SeqLess : less<typename MaxMF<VEC1>::type,typename MaxMF<VEC2>::type> {};


namespace namehider {

using namespace tmptools;


template<typename VA, typename ICW>
struct InnerCheck : fold<VA,false_,or_<numerical_contains<mpl::_2,ICW>,mpl::_1> > {};

template<int RANK, typename VA>
struct Algorithm
  : fold<Ordinals<RANK>,
	 true_,
	 and_<mpl::_1,
	      InnerCheck<VA,mpl::_2>
	      >
	 > 
{};


} // namehider


template<int RANK, typename VA>
struct CheckMeta : namehider::Algorithm<RANK,VA>
{};


} // composite


#endif // TMP_HELPERS_INCLUDED
