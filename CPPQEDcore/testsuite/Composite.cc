// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Composite.h"

#define BOOST_TEST_MODULE CompositeInstantiationAndConstruction test
#include <boost/test/unit_test.hpp>


using namespace std;
using structure::Free;


struct FreeType1 : Free
{
  FreeType1(size_t d) : Free(d) {};
};


struct FreeType2 : Free
{
  FreeType2(size_t d) : Free(d) {};
};


struct FreeType3 : Free
{
  FreeType3(size_t d) : Free(d) {};
};


typedef structure::Interaction<2> I2;
typedef structure::Interaction<3> I3;
typedef structure::Interaction<4> I4;


BOOST_AUTO_TEST_CASE( ONE_CORRECT_AND_ONE_INCORRECT_COMPOSITE )
{

  FreeType1 f0(3), f2(2), f3(5);
  FreeType2 f1(7);
  FreeType3 f4(6);

  I2 i10(I2::Frees(&f1,&f0)), i23(I2::Frees(&f2,&f3)), i43(I2::Frees(&f4,&f3));
  I3 i042(I3::Frees(&f0,&f4,&f2)), i024(I3::Frees(&f0,&f2,&f4));
  I4 i1204(I4::Frees(&f1,&f2,&f0,&f4));

  makeComposite(Act<1,0>(i10),Act<2,3>(i23),Act<4,3>(i43),Act<0,4,2>(i042),Act<1,2,0,4>(i1204));

  try {
    makeComposite(Act<1,0>(i10),Act<2,3>(i23),Act<4,3>(i43),Act<0,2,4>(i042),Act<1,2,0,4>(i1204));
  } catch(const CompositeConsistencyException& cce) {
    cout<<"Exception caught "<<cce.idx<<' '<<cce.i<<endl;
    BOOST_CHECK(cce.idx==2); BOOST_CHECK(cce.i==1);
  }
   

}


namespace incorrect_acts {

// const composite::result_of::make_vector< >::type acts1;

// Does not compile:
// const composite::result_of::make_vector<Act<1,0>,Act<2,0>,Act<3,4>,Act<1,2,0>,Act<2,1,0>,Act<3,3,0,1> >::type acts2;

} // incorrect_acts


namespace a_correct_composite {

typedef composite::result_of::make_vector<Act<1,0>,Act<2,0>,Act<3,4>,Act<1,2,0>,Act<2,1,0>,Act<3,2,0,1> >::type Acts;

} // a_correct_composite

template class Composite<a_correct_composite::Acts>; // explicit instantiation


namespace an_incorrect_composite {

typedef composite::result_of::make_vector<Act<1,0>,Act<2,0>,Act<3,5>,Act<1,2,0>,Act<2,1,0>,Act<3,2,0,1> >::type Acts;

} // an_incorrect_composite

/*
This now won't compile: (COMPOSITE_not_CONSISTENT)
template class Composite<an_incorrect_composite::Acts>; 
*/

