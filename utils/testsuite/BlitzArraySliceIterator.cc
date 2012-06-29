#include "BlitzArraySliceIterator.h"

#include "Randomized.h"

#define USE_BOOST_PROGRESS_TIMER_PROFILING
#include "Profiling.h"

#define DUMMY_TEST_FUNCTION( A ) void A()

#define BOOST_TEST_MODULE BlitzArraySliceIterator test
#include <boost/test/unit_test.hpp>

#include <boost/fusion/sequence/io.hpp>
#include <boost/fusion/sequence/comparison.hpp>
#include <boost/fusion/functional/invocation/invoke.hpp>
#include <boost/fusion/sequence/intrinsic/at_c.hpp>
namespace mpl=boost::mpl;


#include <boost/lambda/lambda.hpp>
namespace bll=boost::lambda;

#include <sstream>

using namespace blitzplusplus;
using namespace basi;
using namespace std;

typedef tmptools::Vector<3,6,1,9,7> V;
const int RANK=11;

typedef details::TransposerMeta<RANK,V>::type TM;

typedef details::IdxTypes<RANK,V>::type Idx;


typedef TTD_IDXTINY(RANK) IdxTiny;
IdxTiny idx(21,2,4,10,11,3,5,9,23,22,7);


typedef TTD_IDXTINY(RANK-mpl::size<V>::value) VecIdxTiny;
VecIdxTiny filteredIdx(details::FilterOut<RANK,V>(idx));


typedef TTD_DARRAY(11) DArray11;


namespace {

using mpl::at_c;
using mpl::equal;

BOOST_STATIC_ASSERT((at_c<TM,0>::type::value==0));
BOOST_STATIC_ASSERT((at_c<TM,1>::type::value==3));
BOOST_STATIC_ASSERT((at_c<TM,2>::type::value==2));
BOOST_STATIC_ASSERT((at_c<TM,3>::type::value==6));
BOOST_STATIC_ASSERT((at_c<TM,4>::type::value==4));
BOOST_STATIC_ASSERT((at_c<TM,5>::type::value==5));
BOOST_STATIC_ASSERT((at_c<TM,6>::type::value==1));
BOOST_STATIC_ASSERT((at_c<TM,7>::type::value==9));
BOOST_STATIC_ASSERT((at_c<TM,8>::type::value==8));
BOOST_STATIC_ASSERT((at_c<TM,9>::type::value==7));
BOOST_STATIC_ASSERT((at_c<TM,10>::type::value==10));

BOOST_STATIC_ASSERT((equal<TM,     mpl::vector_c<int,0,3,2,6,4,5,1,9,8,7,10> >::value));
BOOST_STATIC_ASSERT((equal<TM,tmptools::Vector  <    0,3,2,6,4,5,1,9,8,7,10> >::value));


using blitz::Range;

BOOST_STATIC_ASSERT((equal<Idx,mpl::vector<int,Range,int,Range,int,int,Range,Range,int,Range,int> >::value));


}





BOOST_AUTO_TEST_CASE( FilterOutTest )
{
  using namespace details;

  BOOST_CHECK(all(filteredIdx==VecIdxTiny(21,4,11,3,23,7)));

}


BOOST_AUTO_TEST_CASE( IdxValueTest )
{
  using namespace std;

  const blitz::Range a(blitz::Range::all());

  const Idx v(21,a,4,a,11,3,a,a,23,a,7);

  // The following strange solution is needed because there is no comparison operation for Ranges
  stringstream s1(stringstream::out), s2(stringstream::out);

  struct IndexerBase : details::IndexerBase<RANK,V>
  {
    static const Idx&
    fillIdxValues(const VecIdxTiny& idx) {return details::IndexerBase<RANK,V>::fillIdxValues(idx);}
  };

  s1<<v; s2<<IndexerBase::fillIdxValues(filteredIdx);

  BOOST_CHECK(s1.str()==s2.str());

}


BOOST_AUTO_TEST_CASE( ArraySlicingTest )
{
  typedef mpl::push_front<Idx,DArray11&>::type IdxExt;

  const blitz::Range a(blitz::Range::all());

  DArray11 array11(22,2,5,1,12,4,4,3,24,2,8);

  IdxExt v(array11,21,a,4,a,11,3,a,a,23,a,7);

  // invoke(&DArray11::operator(),v);
  // this will not work because &DArray11::operator() is so heavily overloaded
  
  using boost::fusion::at_c;

  TTD_DARRAY(5) array5(array11(at_c<1>(v),at_c<2>(v),at_c<3>(v),at_c<4>(v),at_c<5>(v),at_c<6>(v),at_c<7>(v),at_c<8>(v),at_c<9>(v),at_c<10>(v),at_c<11>(v)));

  BOOST_CHECK(all(array5.extent()==TTD_EXTTINY(5)(2,1,4,3,2)));

}


BOOST_AUTO_TEST_CASE( ExampleFromManual )
{
  void actWithA(TTD_CARRAY(5)&);

  static const int RANK=11;

  TTD_CARRAY(RANK) psi;

  tmptools::Vector<3,6,1,9,7> v;
  boost::for_each(blitzplusplus::basi::fullRange(psi,v),actWithA);
}

void actWithA(TTD_CARRAY(5)&) {}



// The following test comes from the old version of this file, and is intended to demonstrate the performance as a function of the arity of slices

namespace basi_performance {

const int nRepetition=1;

TTD_CARRAY(11) array1(5,4,5,4,5,4,5,4,5,4,5), array2(array1.shape()), arrayRes(array1.shape()), arrayOrig(array1.shape());

struct Helper
{
  template<typename V> void operator()(V v)
  {
    PROGRESS_TIMER_IN_POINT(cout);
    for (int i=nRepetition; i; --i) cpputils::for_each(fullRange(array1,v),basi::begin(array2,v),bll::_1*=bll::_2); 
    cout<<"Arity "<<11-mpl::size<V>::value<<": ";
    PROGRESS_TIMER_OUT_POINT("");

    BOOST_CHECK(all(array1==arrayRes)); array1=arrayOrig;

    PROGRESS_TIMER_IN_POINT(cout);
    SlicesData<11,V> slicesData(array1);
    for (int i=nRepetition; i; --i) cpputils::for_each(basi_fast::fullRange(array1,slicesData),basi_fast::begin(array2,slicesData),bll::_1*=bll::_2); 
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

TTD_CARRAY(6) array1(3,2,3,2,3,2), array2(array1.shape());


template<int RANK>
void helper(const TTD_CARRAY(RANK)& a1, const TTD_CARRAY(RANK)& a2, const dcomp* dc1, const dcomp* dc2)
{
  cout<<a1.zeroOffset()<<' '<<a1.shape()<<' '<<a1.stride()<<' '<<a1.ordering()<<' '<<a1.data()-dc1<<endl<<a2.zeroOffset()<<' '<<a2.shape()<<' '<<a2.stride()<<' '<<a2.ordering()<<' '<<a2.data()-dc2<<endl;
    BOOST_CHECK(all(a1==a2));
}


struct Helper
{
  
  template<typename V> void operator()(V v)
  {
    SlicesData<6,V> slicesData(array1);

    cout<<endl<<
      "****************\n"<<
      "*** Slice Arity: "<<mpl::size<V>::value<<endl<<
      "****************\n";
    cpputils::for_each(fullRange(array1,v),basi_fast::begin(array2,slicesData),boost::bind(helper<mpl::size<V>::value>,_1,_2,array1.data(),array2.data())); 
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

  /*
  SlicesData<6,tmptools::Vector<3,0,4,1,5> > data(array1);
  basi_fast::Iterator<6,tmptools::Vector<3,0,4,1,5>,true> iter(data,array2,mpl::false_());
  */
}





// BlitzArraySliceIteratorFast assumes: all storage ascending, all bases zero
