#include "TMP_Tools.h"

#include<boost/mpl/size.hpp>
#include<boost/mpl/at.hpp>

#include<boost/static_assert.hpp>

using namespace tmptools;


BOOST_STATIC_ASSERT((
		     Power<3,4>::value==81
		     ));

BOOST_STATIC_ASSERT((
		     Power<2,18>::value==262144
		     ));

// numerical_contains for ...
BOOST_STATIC_ASSERT((
		     numerical_contains_c<
		     RangeMF<10,1>::type,
		     size_t, // ... another type
		     10
		     >::value
		     ));

// 11 should not be in the range 1..10
BOOST_STATIC_ASSERT((
		     !
		     numerical_contains_c<
		     RangeMF<10,1>::type,
		     int,
		     11
		     >::value
		     ));

BOOST_STATIC_ASSERT((
		     numerical_contains_c<
		     Vector<13,1,87,23,198,2834,21>,
		     size_t,
		     21
		     >::value
		     ));

BOOST_STATIC_ASSERT((
		     !
		     numerical_contains_c<
		     Vector<13,1,87,23,198,2834,21>,
		     size_t,
		     22
		     >::value
		     ));


BOOST_STATIC_ASSERT((
		     IsEvenAssert<22>::value==11
		     ));


typedef Vector<87,28,93,1,23,25,97,345,6> V;

BOOST_STATIC_ASSERT((
		     boost::mpl::size<V>::value==9
		     ));


BOOST_STATIC_ASSERT((
		     boost::mpl::at_c<V,2>::type::value==93
		     ));


typedef pair_c<23,42> P_23_42;


BOOST_STATIC_ASSERT((
		     P_23_42::first==23 && P_23_42::second==42
		     ));


void f()
{
  /*
  // Should not compile:
  typedef tmptools::OrdinalMF<-1>::type ordinalError;
  IsEvenAssert<31> error_IsEvenAssert;
  P_23_42::SanityCheck<24,45>();
  P_23_42::SanityCheck<21,40>();
  pair_c<23,23>();
  */

  pair_c<23,23,false>();

  P_23_42::SanityCheck<21,45>();

}
