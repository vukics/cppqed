// -*- C++ -*-

#define ITER BOOST_PP_ITERATION()

#define RETURN_type Composite<typename composite::result_of::make_list<BOOST_PP_ENUM_PARAMS(ITER,T)>::type>

template<BOOST_PP_ENUM_PARAMS(ITER,typename T)> 
const RETURN_type
makeComposite(BOOST_PP_ENUM_BINARY_PARAMS(ITER,const T,& act) )
{
  return RETURN_type(composite::make_list(BOOST_PP_ENUM_PARAMS(ITER,act)));
}

#undef  RETURN_type

#undef  ITER
