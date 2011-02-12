#include "TMP_Tools.h"

#include<boost/mpl/size.hpp>
#include<boost/mpl/at.hpp>

#include<boost/static_assert.hpp>

#include<iostream>

namespace mpl=boost::mpl;

typedef tmptools::RangeMF<2,3>::type Range;

BOOST_STATIC_ASSERT(
		    mpl::size<Range>::value==2
		    );

BOOST_STATIC_ASSERT(
		    (mpl::at_c<Range,0>::type::value==3)
		    );

typedef tmptools::Vector<> emptyVec;
typedef tmptools::Vector<5> Vec5;
typedef tmptools::Vector<3,2> Vec32;
typedef tmptools::Vector<3,2,9,7,6,66,7,5,6,8,7,56,8,55,4,7,87,6> Vec18;

/*
Should not compile:

typedef tmptools::Vector<3,2,9,-1001,6,66,7,5,6,8,7,56,8,55,4,7,87,6> VecError;

typedef VecError::type VVV;
*/

// using namespace boost::mpl;

int main()
{
  using namespace std;

  BOOST_STATIC_ASSERT(
		      mpl::size<emptyVec>::value==0
		      );

  BOOST_STATIC_ASSERT(
		      mpl::size<Vec5>::value==1
		      );

  BOOST_STATIC_ASSERT(
		      (mpl::at_c<Vec5,0>::type::value==5)
		      );

  BOOST_STATIC_ASSERT(
		      mpl::size<Vec32>::value==2
		      );

  BOOST_STATIC_ASSERT(
		      (mpl::at_c<Vec32,0>::type::value==3)
		      );

  BOOST_STATIC_ASSERT(
		      (mpl::at_c<Vec32,1>::type::value==2)
		      );

  BOOST_STATIC_ASSERT(
		      mpl::size<Vec18>::value==18
		      );

  BOOST_STATIC_ASSERT(
		      (mpl::at_c<Vec18,13>::type::value==55)
		      );

  BOOST_STATIC_ASSERT(
		      mpl::bool_< (1>0) >::value
		      );
  BOOST_STATIC_ASSERT(
		      (boost::mpl::if_c<(19 < 0),mpl_::void_,mpl_::int_<19> >::type::value==19)
		      );
}
