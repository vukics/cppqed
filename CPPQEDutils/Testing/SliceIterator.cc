// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "BlitzArray.h"
#include "SliceIterator.tcc"

// #include "Random.h"

#define BOOST_TEST_MODULE BlitzArraySliceIterator test
#include <boost/test/unit_test.hpp>

#include <boost/fusion/sequence/io.hpp>
#include <boost/fusion/sequence/comparison.hpp>
#include <boost/fusion/functional/invocation/invoke.hpp>
#include <boost/fusion/sequence/intrinsic/at_c.hpp>
namespace mpl=boost::mpl;

#include <sstream>

using namespace cppqedutils::sliceiterator;
using namespace std;

using V=tmptools::Vector<3,6,1,9,7>;

const int RANK=11;

using TM=details::TransposerMeta_t<RANK,V>;

using FullIdx=details::IndexerBase<RANK,V>::Idx; // This contains Range types at the appropriate locations


const IdxTiny<RANK> idx(21,2,4,10,11,3,5,9,23,22,7);


const VecIdxTiny<RANK,V> filteredIdx(filterOut<RANK,V>(idx));

using DArray11=DArray<11>;


namespace {

using mpl::at_c, mpl::equal;

static_assert(at_c<TM,0>::type::value==0);
static_assert(at_c<TM,1>::type::value==3);
static_assert(at_c<TM,2>::type::value==2);
static_assert(at_c<TM,3>::type::value==6);
static_assert(at_c<TM,4>::type::value==4);
static_assert(at_c<TM,5>::type::value==5);
static_assert(at_c<TM,6>::type::value==1);
static_assert(at_c<TM,7>::type::value==9);
static_assert(at_c<TM,8>::type::value==8);
static_assert(at_c<TM,9>::type::value==7);
static_assert(at_c<TM,10>::type::value==10);

static_assert(equal<TM,     mpl::vector_c<int,0,3,2,6,4,5,1,9,8,7,10> >::value);
static_assert(equal<TM,tmptools::Vector  <    0,3,2,6,4,5,1,9,8,7,10> >::value);


using blitz::Range;

static_assert(equal<FullIdx,mpl::vector<int,Range,int,Range,int,int,Range,Range,int,Range,int> >::value);


}


BOOST_AUTO_TEST_CASE( FilterOutTest )
{
  BOOST_CHECK(all(filteredIdx==VecIdxTiny<RANK,V>(21,4,11,3,23,7)));
}


BOOST_AUTO_TEST_CASE( IdxValueTest )
{
  using namespace std;

  const blitz::Range a(blitz::Range::all());

  const FullIdx v(21,a,4,a,11,3,a,a,23,a,7);

  // The following strange solution is needed because there is no comparison operation for Ranges
  stringstream s1(stringstream::out), s2(stringstream::out);

  s1<<v; s2<<details::IndexerBase<RANK,V>::fillIdxValues(filteredIdx);

  BOOST_CHECK(s1.str()==s2.str());

}


BOOST_AUTO_TEST_CASE( ArraySlicingTest )
{
  const auto a{blitz::Range::all()};

  DArray11 array11{22,2,5,1,12,4,4,3,24,2,8};

  FullIdx v{21,a,4,a,11,3,a,a,23,a,7};

  using boost::fusion::at_c;

  BOOST_CHECK(
    all(DArray<5>{array11(at_c<0>(v),at_c<1>(v),at_c<2>(v),at_c<3>(v),at_c<4>(v),at_c<5>(v),at_c<6>(v),at_c<7>(v),at_c<8>(v),at_c<9>(v),at_c<10>(v))}.extent()
        ==
        ExtTiny<5>(2,1,4,3,2)));

}


BOOST_AUTO_TEST_CASE( ExampleFromManual )
{
  void actWithA(CArray<5>&);

  CArray<RANK> psi;

  for (auto p : fullRange<V>(psi) ) actWithA(p);
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
