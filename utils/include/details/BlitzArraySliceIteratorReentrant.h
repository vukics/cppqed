// -*- C++ -*-

template<typename A, typename V>
inline
const RETURN_type1(true )
begin(ADD_PARAMETER ADD_parameter const A& array, V)
{return RETURN_type1(true )(ADD_parameter array,boost::mpl::false_());}

template<typename A, typename V>
inline
const RETURN_type1(true )
end  (ADD_PARAMETER ADD_parameter const A& array, V) 
{return RETURN_type1(true )(ADD_parameter array,boost::mpl:: true_());}

template<typename A, typename V>
inline
const RETURN_type1(false)
begin(ADD_PARAMETER ADD_parameter       A& array, V) 
{return RETURN_type1(false)(ADD_parameter array,boost::mpl::false_());}

template<typename A, typename V>
inline
const RETURN_type1(false)
end  (ADD_PARAMETER ADD_parameter       A& array, V) 
{return RETURN_type1(false)(ADD_parameter array,boost::mpl:: true_());}


#define RETURN_type2(CONST) boost::iterator_range<RETURN_type1(CONST)>

template<typename A, typename V>
inline
const RETURN_type2(true ) 
fullRange(ADD_PARAMETER ADD_parameter const A& array, V v) 
{return RETURN_type2(true )(blitzplusplus::NS_NAME::begin<A,V>(ADD_parameter array,v),blitzplusplus::NS_NAME::end<A,V>(ADD_parameter array,v));}

template<typename A, typename V>
inline
const RETURN_type2(false) 
fullRange(ADD_PARAMETER ADD_parameter       A& array, V v) 
{return RETURN_type2(false)(blitzplusplus::NS_NAME::begin<A,V>(ADD_parameter array,v),blitzplusplus::NS_NAME::end<A,V>(ADD_parameter array,v));}

#undef  RETURN_type2
#undef  RETURN_type1
#undef  ADD_PARAMETER
#undef  ADD_parameter
#undef  NS_NAME
