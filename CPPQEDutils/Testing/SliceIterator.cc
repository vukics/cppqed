// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

#include "SliceIterator.h"

#include "BlitzArray.h"

#include "TMP_Tools.h"

// #include "Random.h"

#define BOOST_TEST_MODULE BlitzArraySliceIterator test
#include <boost/test/unit_test.hpp>

#include <sstream>

using namespace cppqedutils::sliceiterator;
using namespace hana::literals;
using namespace std;

using RetainedAxes=tmptools::Vector<3,6,1,9,7>;

const int RANK=11;


constexpr auto ti = transposingIndeces<RANK,RetainedAxes>();

const auto fullIdx=slicingTuple<RANK,RetainedAxes>(); // This contains Ranges at the appropriate locations


const IdxTiny<RANK> idx(21,2,4,10,11,3,5,9,23,22,7);


const auto filteredIdx(filterOut<RANK,RetainedAxes>(idx));

using DArray11=DArray<11>;


BOOST_HANA_CONSTANT_ASSERT( ti[0_c] == 0_c );
BOOST_HANA_CONSTANT_ASSERT( ti[1_c] == 3_c );
BOOST_HANA_CONSTANT_ASSERT( ti[2_c] == 2_c );
BOOST_HANA_CONSTANT_ASSERT( ti[3_c] == 6_c );
BOOST_HANA_CONSTANT_ASSERT( ti[4_c] == 4_c );
BOOST_HANA_CONSTANT_ASSERT( ti[5_c] == 5_c );
BOOST_HANA_CONSTANT_ASSERT( ti[6_c] == 1_c );
BOOST_HANA_CONSTANT_ASSERT( ti[7_c] == 9_c );
BOOST_HANA_CONSTANT_ASSERT( ti[8_c] == 8_c );
BOOST_HANA_CONSTANT_ASSERT( ti[9_c] == 7_c );
BOOST_HANA_CONSTANT_ASSERT( ti[10_c] == 10_c );

static_assert( ti[9_c] == 7_c );

static_assert( ti == tmptools::vector<0,3,2,6,4,5,1,9,8,7,10> );

namespace {

using blitz::Range;

static_assert( hana::type_c<std::decay_t<decltype(fullIdx[2_c])>> == hana::type_c<int> );
static_assert( hana::type_c<std::decay_t<decltype(fullIdx[6_c])>> == hana::type_c<Range> );

static_assert( std::is_same_v<decltype(fullIdx),const hana::tuple<int,Range,int,Range,int,int,Range,Range,int,Range,int> > );

}


BOOST_AUTO_TEST_CASE( FilterOutTest )
{
  BOOST_CHECK(all(filteredIdx==decltype(filteredIdx){21,4,11,3,23,7}));
}


BOOST_AUTO_TEST_CASE( TransposeTest )
{
  DArray11 array11{22,2,5,1,12,4,4,3,24,2,8};
  transpose<RetainedAxes>(array11);
  std::cout<<array11.extent();
  BOOST_CHECK(all( array11.extent() == ExtTiny<11>{22,1,5,4,12,4,2,2,24,3,8} ));
}


/*
BOOST_AUTO_TEST_CASE( IdxValueTest )
{
  using namespace std;

  const blitz::Range a(blitz::Range::all());

  const FullIdx v(21,a,4,a,11,3,a,a,23,a,7);

  // The following strange solution is needed because there is no comparison operation for Ranges
  stringstream s1(stringstream::out), s2(stringstream::out);

  s1<<v; s2<<details::IndexerBase<RANK,RetainedAxes>::fillIdxValues(filteredIdx);

  BOOST_CHECK(s1.str()==s2.str());

}
*/


BOOST_AUTO_TEST_CASE( ArraySlicingTest )
{
  const auto a{blitz::Range::all()};

  DArray11 array11{22,2,5,1,12,4,4,3,24,2,8};

  {
    using namespace hana::literals;
    
    decltype(fullIdx) v{21,a,4,a,11,3,a,a,23,a,7};    
    
    BOOST_CHECK( all( DArray<5>{array11(v[0_c],v[1_c],v[2_c],v[3_c],v[4_c],v[5_c],v[6_c],v[7_c],v[8_c],v[9_c],v[10_c])}.extent() == ExtTiny<5>{2,1,4,3,2} ) );
  
  }
  
  DArray<5> array5;
  
  {
    Slicer<11,RetainedAxes>::_<DArray>(array11,array5,IdxTiny<6>{18,3,8,1,15,5});
    BOOST_CHECK( all( array5.extent() == ExtTiny<5>{2,1,4,3,2} ) );
  }
  
}


BOOST_AUTO_TEST_CASE( ExampleFromManual )
{
  void actWithA(CArray<5>&);

  CArray<RANK> psi;

  for (auto p : fullRange<RetainedAxes>(psi) ) actWithA(p);
}

void actWithA(CArray<5>&) {}


// The following test comes from the old version of this file, and is intended to demonstrate the performance as a function of the arity of slices
/*
namespace basi_performance {

const int nRepetition=1;

CArray<11> array1(5,4,5,4,5,4,5,4,5,4,5), array2(array1.shape()), arrayRes(array1.shape()), arrayOrig(array1.shape());

struct Helper
{
  template<typename V> void operator()(V)
  {
    PROGRESS_TIMER_IN_POINT(cout);
    for (int i=nRepetition; i; --i) cppqedutils::for_each(fullRange<V>(array1),basi::begin<V>(array2),bll::_1*=bll::_2); 
    cout<<"Arity "<<11-mpl::size<V>::value<<": ";
    PROGRESS_TIMER_OUT_POINT("");

    BOOST_CHECK(all(array1==arrayRes)); array1=arrayOrig;

    PROGRESS_TIMER_IN_POINT(cout);
    SlicesData<11,V> slicesData(array1);
    for (int i=nRepetition; i; --i) cppqedutils::for_each(basi_fast::fullRange(array1,slicesData),basi_fast::begin(array2,slicesData),bll::_1*=bll::_2); 
    cout<<"Fast. Arity "<<11-mpl::size<V>::value;
    PROGRESS_TIMER_OUT_POINT("");

    BOOST_CHECK(all(array1==arrayRes)); array1=arrayOrig;
  }

};

} // basi_performance


BOOST_AUTO_TEST_CASE( BASI_Performance )
{

  using namespace randomized; using namespace basi_performance;

  fillWithRandom(array2,fillWithRandom(array1));

  arrayRes=arrayOrig=array1;

  PROGRESS_TIMER_IN_POINT(cout);
  for (int i=nRepetition; i; --i) arrayRes*=array2;
  PROGRESS_TIMER_OUT_POINT("\nBlitz internal multiplication");

  mpl::for_each<mpl::vector<
    tmptools::Vector<9,3,6,0,10,4,7,8,1,5,2>,
    tmptools::Vector<9,3,6,0,10,4,7,8,1,5>,
    tmptools::Vector<9,3,6,0,10,4,7,8,1>,
    tmptools::Vector<9,3,6,0,10,4,7,8>,
    tmptools::Vector<9,3,6,0,10,4,7>,
    tmptools::Vector<9,3,6,0,10,4>,
    tmptools::Vector<9,3,6,0,10>,
    tmptools::Vector<9,3,6,0>,
    tmptools::Vector<9,3,6>,
    tmptools::Vector<9,3>,
    tmptools::Vector<9>
      > >(Helper());

}





namespace basi_monitor {

CArray<6> array1(3,2,3,2,3,2), array2(array1.shape());


template<int RANK>
void helper(const CArray<RANK>& a1, const CArray<RANK>& a2, const dcomp* dc1, const dcomp* dc2)
{
  cout<<a1.zeroOffset()<<' '<<a1.shape()<<' '<<a1.stride()<<' '<<a1.ordering()<<' '<<a1.data()-dc1<<endl<<a2.zeroOffset()<<' '<<a2.shape()<<' '<<a2.stride()<<' '<<a2.ordering()<<' '<<a2.data()-dc2<<endl;
  BOOST_CHECK(all(a1==a2));
}


struct Helper
{
  
  template<typename V> void operator()(V)
  {
    SlicesData<6,V> slicesData(array1);

    cout<<endl<<
      "****************\n"<<
      "*** Slice Arity: "<<mpl::size<V>::value<<endl<<
      "****************\n";
    cppqedutils::for_each(fullRange<V>(array1),basi_fast::begin(array2,slicesData),boost::bind(helper<mpl::size<V>::value>,_1,_2,array1.data(),array2.data())); 
  }

};

} // basi_monitor



BOOST_AUTO_TEST_CASE( BASI_Monitor )
{

  using namespace basi_monitor;

  randomized::fillWithRandom(array1);

  array2=array1;

  mpl::for_each<mpl::vector<
  tmptools::Vector<3,0,4,1,5>,
  tmptools::Vector<3,0,4,1>,
  tmptools::Vector<3,0,4>,
  tmptools::Vector<3,0>,
  tmptools::Vector<3>
    > > (Helper());

  
  // SlicesData<6,tmptools::Vector<3,0,4,1,5> > data(array1);
  // basi_fast::Iterator<6,tmptools::Vector<3,0,4,1,5>,true> iter(data,array2,mpl::false_());

}
*/


// BlitzArraySliceIteratorFast assumes: all storage ascending, all bases zero
