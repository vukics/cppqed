// -*- C++ -*-

#define ADD_argument BOOST_PP_IIF(IS_SD,sd,) BOOST_PP_COMMA_IF(IS_SD)

template<typename A, typename SD_V>
inline
const RETURN_type1(true )
begin(const A& array, BOOST_PP_IIF(IS_SD,const SD_V&,SD_V) BOOST_PP_IIF(IS_SD,sd,) )
{return RETURN_type1(true )(array, ADD_argument boost::mpl::false_());}

template<typename A, typename SD_V>
inline
const RETURN_type1(true )
end  (const A& array, BOOST_PP_IIF(IS_SD,const SD_V&,SD_V) BOOST_PP_IIF(IS_SD,sd,) )
{return RETURN_type1(true )(array, ADD_argument boost::mpl:: true_());}

template<typename A, typename SD_V>
inline
const RETURN_type1(false)
begin(      A& array, BOOST_PP_IIF(IS_SD,const SD_V&,SD_V) BOOST_PP_IIF(IS_SD,sd,) )
{return RETURN_type1(false)(array, ADD_argument boost::mpl::false_());}

template<typename A, typename SD_V>
inline
const RETURN_type1(false)
end  (      A& array, BOOST_PP_IIF(IS_SD,const SD_V&,SD_V) BOOST_PP_IIF(IS_SD,sd,) ) 
{return RETURN_type1(false)(array, ADD_argument boost::mpl:: true_());}


#define RETURN_type2(CONST) boost::iterator_range<RETURN_type1(CONST)>

template<typename A, typename SD_V>
inline
const RETURN_type2(true ) 
fullRange(const A& array, BOOST_PP_IIF(IS_SD,const SD_V&,SD_V) sd_v ) 
{return RETURN_type2(true )(blitzplusplus::NS_NAME::begin<A,SD_V>(array,sd_v),blitzplusplus::NS_NAME::end<A,SD_V>(array,sd_v));}

template<typename A, typename SD_V>
inline
const RETURN_type2(false) 
fullRange(      A& array, BOOST_PP_IIF(IS_SD,const SD_V&,SD_V) sd_v ) 
{return RETURN_type2(false)(blitzplusplus::NS_NAME::begin<A,SD_V>(array,sd_v),blitzplusplus::NS_NAME::end<A,SD_V>(array,sd_v));}


#undef  RETURN_type2
#undef  ADD_argument
#undef  IS_SD
#undef  RETURN_type1
#undef  NS_NAME
