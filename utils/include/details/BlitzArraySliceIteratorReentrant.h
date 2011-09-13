// -*- C++ -*-

template<typename A, typename V_S>
inline
const RETURN_type1(true )
begin(const A& array, V_S )
{return RETURN_type1(true )(array, boost::mpl::false_());}

template<typename A, typename V_S>
inline
const RETURN_type1(true )
end  (const A& array, V_S )
{return RETURN_type1(true )(array, boost::mpl:: true_());}

template<typename A, typename V_S>
inline
const RETURN_type1(false)
begin(      A& array, V_S )
{return RETURN_type1(false)(array, boost::mpl::false_());}

template<typename A, typename V_S>
inline
const RETURN_type1(false)
end  (      A& array, V_S ) 
{return RETURN_type1(false)(array, boost::mpl:: true_());}


#define RETURN_type2(CONST) boost::iterator_range<RETURN_type1(CONST)>

template<typename A, typename V_S>
inline
const RETURN_type2(true ) 
fullRange(const A& array, V_S v_s ) 
{return RETURN_type2(true )(blitzplusplus::NS_NAME::begin<A,V_S>(array,v_s),blitzplusplus::NS_NAME::end<A,V_S>(array,v_s));}

template<typename A, typename V_S>
inline
const RETURN_type2(false) 
fullRange(      A& array, V_S v_s ) 
{return RETURN_type2(false)(blitzplusplus::NS_NAME::begin<A,V_S>(array,v_s),blitzplusplus::NS_NAME::end<A,V_S>(array,v_s));}


#undef  RETURN_type2
#undef  RETURN_type1
#undef  NS_NAME
