// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/**
 * TODO: Consider whether unit tests could be included into the source files directly.
 * This is probably best done in .cc files, of which we do not have too many left …
 */
#include "MultiArrayComplex.h"
#include "Random.h"

#define BOOST_TEST_MODULE MultiArray test
#include <boost/test/unit_test.hpp>

#include <iostream>

using namespace cppqedutils;

constexpr auto ra20741 = retainedAxes<2,0,7,4,1>;

using CArray9=MultiArray<dcomp,9>;

CArray9 array{{5,2,3,2,4,3,2,4,6}, [] (size_t e) {
  CArray9::StorageType res(e);
  std::uniform_real_distribution d{-1.0,1.0};
  std::mt19937_64 randomEngine{1001};
  std::ranges::generate(res, [&]() {
    double re{d(randomEngine)}, im{d(randomEngine)};
    return std::complex{re,im};});
  return res;
}};

Extents<5> filteredExtents{multiarray::filterIn<ra20741>(array.extents)}/*,
  filteredStrides{multiarray::filterIn<ra20741>(array.strides)}*/;


BOOST_AUTO_TEST_CASE( MULTIARRAY_EXTENTS_INDEXING , * boost::unit_test::tolerance(1e-6) )
{
  dcomp element{array(1,1,2,1,0,1,1,3,2)};
  std::cerr<<array.dataView.size()<<std::endl<<element<<std::endl<<json( filteredExtents )<<std::endl;
  BOOST_TEST ( array.dataView.size() == 34560 ) ;
  BOOST_TEST ( std::abs(element-dcomp{-0.864561,0.667605}) == 0 ) ; // direct complex comparison doesn’t work for some reason :-O
  BOOST_TEST ( ( filteredExtents == Extents<5>{3,5,4,4,2} ) ) ;
}


auto sliceRangeTester(const auto& range)
{
  auto iter{range.begin()};
  auto res0=iter->extents;

  for (size_t i=0; i<10; (++i, ++iter));

  auto res1=(*iter)(2,2,1,3,0);

  auto rangeRecurse{sliceRange<retainedAxes<1,3>>(*iter)};

  auto res2=rangeRecurse.begin()->extents;

  auto iter2{rangeRecurse.begin()};
  for (size_t i=0; i<3; (++i, ++iter2));

  auto res3=(*iter2)(3,2);

  std::cerr<<json(res0)<<std::endl<<res1<<std::endl<<json(res2)<<std::endl<<res3<<std::endl;

  return std::make_tuple(res0,res1,res2,res3);
}


void testEvaluator(const auto& res)
{
  BOOST_TEST ( ( std::get<0>(res) == Extents<5>{3,5,4,4,2} ) ) ;
  BOOST_TEST (  std::abs(std::get<1>(res)-dcomp(0.424108,-0.815181)) == 0 ) ;
  BOOST_TEST ( ( std::get<2>(res) == Extents<2>{5,4} ) ) ;
  BOOST_TEST (  std::abs(std::get<3>(res)-dcomp(0.608949,0.275701)) == 0 ) ;
}


BOOST_AUTO_TEST_CASE( MULTIARRAY_SLICING_SIMPLE_RANGE , * boost::unit_test::tolerance(1e-6) )
{
  testEvaluator(sliceRangeTester(sliceRangeSimple<ra20741>(array)));
}


BOOST_AUTO_TEST_CASE( MULTIARRAY_SLICING_PROXY_RANGE_OWNING , * boost::unit_test::tolerance(1e-6) )
{
  testEvaluator(sliceRangeTester(sliceRange<ra20741>(array)));
}


BOOST_AUTO_TEST_CASE( MULTIARRAY_SLICING_PROXY_RANGE_REFERENCING , * boost::unit_test::tolerance(1e-6) )
{
  auto offsets{calculateSlicesOffsets<ra20741>(array.extents)};
  testEvaluator(sliceRangeTester(sliceRange<ra20741>(array,offsets)));
}


/// TODO: test directProduct
/*
  {
    MultiArray<dcomp,4> a1{{5,2,3,2}};
    MultiArray<dcomp,5> a2{{4,3,2,4,6}};

    {
      std::uniform_real_distribution d{-1.0,1.0};
      std::mt19937_64 re{1001};
      std::ranges::generate(a1.dataStorage(), [&]() {return std::complex{d(re),d(re)};} );
      std::ranges::generate(a2.dataStorage(), [&]() {return std::complex{d(re),d(re)};} );
    }

    auto res{directProduct<4,5>(a1,a2)};

  }

  {
    const CArray9& carray=array;

    auto v{vectorize(array)};
    auto cv{vectorize(carray)};
  }

  {
    MultiArray<dcomp,8> array1{{5,2,3,2,5,2,3,2}};
    const MultiArray<dcomp,8>& carray1 = array1;
    auto m{matricize(array1)};
    auto cm{matricize(carray1)};
  }

}

*/
