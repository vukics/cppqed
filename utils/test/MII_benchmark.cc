#include "impl/BlitzArraySliceIterator.tcc"
#include "Randomized.h"

#include <boost/progress.hpp>

#include <list>

using namespace randomized;
using namespace std;

typedef TTD_IDXTINY(10) Size    ;
typedef TTD_IDXTINY( 6) DummyIdx;
typedef cpputils::MultiIndexIterator<6> MII6;
typedef tmptools::Vector<6,2,5,7> SliceVec;
typedef blitzplusplus::basi::Indexer<10,SliceVec> Indexer;
typedef blitzplusplus::basi::Iterator<10,SliceVec,true> BASI;
typedef blitzplusplus::SlicesData<10,SliceVec> SlicesData;
typedef blitzplusplus::basi_fast::Iterator<10,SliceVec,true> BASI_FAST;

const Size size(3,4,2,5,3,5,6,4,2,5);

const SliceVec sliceVec=SliceVec();

const size_t nRepeat=
#ifndef NDEBUG
  1
#else
  10000
#endif
  ;

int main()
{
  TTD_CARRAY(10) array(size);
  TTD_CARRAY( 4) arrayRes/*(6,2,5,4)*/;

  fillWithRandom(array);

  { // 1

  list<DummyIdx> dummyIdcs;

  const MII6
    begin(DummyIdx(0,0,0,0,0,0),DummyIdx(2,3,4,2,1,4),boost::mpl::false_()),
    end  (DummyIdx(0,0,0,0,0,0),DummyIdx(2,3,4,2,1,4),boost::mpl:: true_());

  copy(begin,end,back_inserter(dummyIdcs));

  // cout<<*begin<<endl<<*end<<endl;

  dcomp res;

  {
    boost::progress_timer t;
    for (size_t count=0; count<nRepeat; ++count) {
      for (MII6 i=begin; i!=end; ++i)
	res=Indexer::index(array,arrayRes,*i)(1,2,4,3);
    }
  }

  } // 1

  { // 2

  dcomp res;

  const BASI
    begin(array,boost::mpl::false_()),
    end  (array,boost::mpl:: true_());

  {
    boost::progress_timer t;
    for (size_t count=0; count<nRepeat; ++count) {
      for (BASI i=begin; i!=end; ++i)
	res=(*i)(2,1,4,3);
    }
  }

  } // 2

  { // 3

  const SlicesData sd(array);

  const BASI_FAST
    begin(array,sd,boost::mpl::false_()),
    end  (array,sd,boost::mpl:: true_());

  {
    dcomp res;

    boost::progress_timer t;
    for (size_t count=0; count<nRepeat; ++count) {
      for (BASI_FAST i=begin; i!=end; ++i)
	res=(*i)(2,1,4,3);
    }
  }

  } // 3

}
