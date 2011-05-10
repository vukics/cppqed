namespace blitzplusplus {

namespace basi {

namespace namehider {

using namespace boost::mpl;
namespace mpl=boost::mpl;
using namespace tmptools;

template<int RANK, typename V>
struct Algorithm 
  : fold<typename OrdinalMF<RANK>::type,
	 pair<vector_c<int>,typename boost::mpl::begin<V>::type>,
	 pair<push_back<first<mpl::_1>,
			if_<numerical_contains<V,mpl::_2>,
			    deref<second<mpl::_1> >,
			    mpl::_2
			    >
			>,
	      if_<numerical_contains<V,mpl::_2>,
		  next<second<mpl::_1> >,
		  second<mpl::_1>
		  >
	      >
	 >
{};

} // namehider


template<int RANK, typename V>
struct TransposerMeta : boost::mpl::first<typename namehider::Algorithm<RANK,V>::type>
{};

} // basi

} // blitzplusplus