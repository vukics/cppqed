// -*- C++ -*-

#define rank BOOST_PP_ITERATION()

#define INDEXER_print1(z,m,data) boost::mpl::at_c<TM,m>::type::value
#define INDEXER_print2(z,m,data) boost::fusion::at_c<m>(Base::cache_)

namespace blitzplusplus {
namespace basi {

template<typename V> struct Transposer<rank,V>
{
  typedef typename details::TransposerMeta<rank,V>::type TM;

  typedef TTD_CARRAY(rank) Array;

  static Array& transpose(Array& array)
  {
    array.transposeSelf(BOOST_PP_ENUM(rank,INDEXER_print1,~) );
    return array;
  }

};


template<typename V> struct Indexer<rank,V> : Transposer<rank,V>, private details::IndexerBase<rank,V>
{
  typedef details::IndexerBase<rank,V> Base;

  typedef TTD_CARRAY(rank)  Array   ;
  typedef TTD_RES_CARRAY(V) ArrayRes;

  static ArrayRes& index(Array& array, ArrayRes& arrayRes, const typename Base::VecIdxTiny& idx)
  {
    Base::fillIdxValues(idx);
    arrayRes.reference(array(BOOST_PP_ENUM(rank,INDEXER_print2,~) ) );
    return arrayRes;
  }

};

}
}

#undef INDEXER_print1
#undef INDEXER_print2

#undef rank
