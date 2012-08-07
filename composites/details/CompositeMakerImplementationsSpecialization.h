// -*- C++ -*-

#define ITER BOOST_PP_ITERATION()

namespace composite {


namespace result_of {


template<BOOST_PP_ENUM_PARAMS(ITER,typename A)>
struct Make<BOOST_PP_ENUM_PARAMS(ITER,A) BOOST_PP_ENUM_TRAILING(BOOST_PP_SUB(FUSION_MAX_VECTOR_SIZE,ITER),DEFAULT_print,~) >
  : boost::mpl::identity<Composite<typename make_list<BOOST_PP_ENUM_PARAMS(ITER,A) >::type> >
{};


} // result_of


#define RETURN_type typename result_of::Make<BOOST_PP_ENUM_PARAMS(ITER,A) >::type

template<BOOST_PP_ENUM_PARAMS(ITER,typename A)> 
const RETURN_type
make(BOOST_PP_ENUM_BINARY_PARAMS(ITER,const A,& act) )
{
  return RETURN_type(make_list(BOOST_PP_ENUM_PARAMS(ITER,act)));
}

#undef  RETURN_type


} // composite


#undef  ITER
