// -*- C++ -*-
#ifndef _TRANSFORMATION_IMPL_H
#define _TRANSFORMATION_IMPL_H

#include <boost/fusion/adapted/mpl.hpp>
#include <boost/fusion/tuple.hpp>
#include <boost/fusion/view/zip_view.hpp>
#include <boost/fusion/view/filter_view.hpp>
#include <boost/mpl/count_if.hpp>



namespace quantumdata {

namespace transformation {


namespace specializations {

using namespace linalg;

void performTransformation(const CMatrix      & trafo, const CVector      & in, CVector      & out);

void performTransformation(const TTD_CARRAY(4)& trafo, const TTD_CARRAY(2)& in, TTD_CARRAY(2)& out);

// ... Et caetera.

} // specializations



template<int TWO_TIMES_RANK>
void
Traits<TTD_CARRAY(TWO_TIMES_RANK)>::transform(const TTD_CARRAY(TWO_TIMES_RANK)& trafo, const StateVectorLow& in, StateVectorLow& out)
{
  specializations::performTransformation(trafo,in,out);
}



///////////////////////////
// Composite transformation
///////////////////////////


namespace namehider {


using namespace boost::mpl;
namespace mpl=boost::mpl;


template<typename TRAFO >
struct RankThroughTraits :  int_< Traits<TRAFO>::N_RANK     >
{};


template<typename TRAFO >
struct NotIdentity : true_
{};

template<int RANK>
struct NotIdentity<Identity<RANK> > : false_
{};


template<typename TRAFOS>
struct Algorithm 
  : fold<TRAFOS,
	 pair<tmptools::Vector<>,int_<0> >,
	 pair<if_<NotIdentity<mpl::_2>,
		  push_back<first<mpl::_1>,second<mpl::_1> >,
		  first<mpl::_1>
		  >,
	      plus<second<mpl::_1>,
		   RankThroughTraits<mpl::_2>
		   >
	      >
	 >
{};





} // namehider




template<typename TRAFOS>
class Composite
{
public:
  static const int N_RANK=boost::mpl::fold<TRAFOS,
					   boost::mpl::int_<0>,
					   boost::mpl::plus<boost::mpl::_1,
							    namehider::RankThroughTraits<boost::mpl::_2>
							    >
					   >::type::value;

  typedef TRAFOS TrafoTypes;
  
private:
  typedef boost::fusion::filter_view<const TrafoTypes,namehider::NotIdentity<boost::mpl::_> > TrafoTypesFilteredView;

public:
  typedef typename boost::fusion::result_of::as_vector<TrafoTypesFilteredView>::type FilteredTrafoTypes;
  // Without Identities. It is important to convert back to a vector.

  typedef typename Types<N_RANK>::StateVectorLow StateVectorLow;

  Composite(const TrafoTypes& t): trafos_(t), filteredTrafos_(TrafoTypesFilteredView(trafos_)) {};

  const TrafoTypes& getTrafos() const {return trafos_;}


  void transform(const StateVectorLow& in, StateVectorLow& out) const
  {
    doTransform(in,out,boost::mpl::int_<boost::mpl::count_if<TRAFOS,namehider::NotIdentity<boost::mpl::_> >::value>());
  }


private:

  const TrafoTypes trafos_;

  // implementation helpers:

  const FilteredTrafoTypes filteredTrafos_;

  DEFINE_TYPED_STATIC_CONST( typename boost::mpl::first<typename namehider::Algorithm<TRAFOS>::type>::type , AuxiliaryVector , auxiliaryVector_ ) ;

  const boost::fusion::zip_view<boost::fusion::tuple<const FilteredTrafoTypes&, const AuxiliaryVector&> >
  zip() const 
  {
    return tie(filteredTrafos_,auxiliaryVector_); // implicit conversion to return type
  }

  template<int COUNT>
  void doTransform(const StateVectorLow& in, StateVectorLow& out, boost::mpl::int_<COUNT>) const
  {
    using namespace boost::fusion;
    StateVectorLow buf(out.shape());
    boost::fusion::fold(zip(),in,TransformSubsystem(out,buf));
  }

  void doTransform(const StateVectorLow& in, StateVectorLow& out, boost::mpl::int_<  1  >) const
  {
    using namespace boost::fusion;
    boost::fusion::fold(zip(),in,TransformSubsystem(out,out));
  }

  void doTransform(const StateVectorLow& in, StateVectorLow& out, boost::mpl::int_<  0  >) const
  {
    out=in;
  }


  class TransformSubsystem
  {
  public:
    typedef const StateVectorLow& result_type;

    TransformSubsystem(StateVectorLow& out, StateVectorLow& buf)
        : out_(out), buf_(buf) {}


    template<typename TUPLE>
    const StateVectorLow& doIt(const StateVectorLow& in, TUPLE, boost::mpl:: true_)
    {
      return in;
    }

#define TYPEDEF_TRAFO typedef typename tmptools::RemoveConstReference<typename boost::fusion::result_of::at_c<TUPLE,0>::type>::type TRAFO;

    template<typename TUPLE>
    const StateVectorLow& doIt(const StateVectorLow& in, TUPLE tuple, boost::mpl::false_)
    {
      using namespace boost::fusion; using namespace blitzplusplus; using basi::fullRange;

      TYPEDEF_TRAFO ;

      typedef typename tmptools::RangeMF<Traits<TRAFO>::N_RANK,boost::remove_reference<typename result_of::at_c<TUPLE,1>::type>::type::value>::type Range;

      cpputils::for_each(fullRange(in,Range()),basi::begin(buf_,Range()),bind(&Traits<TRAFO>::transform,at_c<0>(tuple),_1,_2));

      blitz::swap(buf_,out_); return out_;
    }

    template<typename TUPLE>
    const StateVectorLow& operator()(const StateVectorLow& in, TUPLE tuple)
    {
      TYPEDEF_TRAFO ;

      return doIt(in,tuple,boost::mpl::bool_<!namehider::NotIdentity<TRAFO>::value>());
    }

#undef TYPEDEF_TRAFO

    template<typename TUPLE>
    const StateVectorLow& operator()(TUPLE tuple, const StateVectorLow& in)
    {
      return operator()(in,tuple);
    }


  private:
    StateVectorLow& out_, buf_;

  }; // TransformSubsystems


}; // Composite



template<typename TRAFOS>
const typename Composite<TRAFOS>::AuxiliaryVector Composite<TRAFOS>::auxiliaryVector_=Composite<TRAFOS>::AuxiliaryVector();


#define COMPOSITE_TRAFO Composite<typename boost::fusion::result_of::as_vector<typename boost::fusion::result_of::join<const typename Traits<TRAFO1>::TrafoTypes,const typename Traits<TRAFO2>::TrafoTypes>::type>::type>

template<typename TRAFO1, typename TRAFO2>
struct Compose : boost::mpl::identity<COMPOSITE_TRAFO>
{
  static const COMPOSITE_TRAFO compose(const TRAFO1& t1, const TRAFO2& t2)
  {
    return COMPOSITE_TRAFO(boost::fusion::as_vector(boost::fusion::join(Traits<TRAFO1>::getTrafos(t1),Traits<TRAFO2>::getTrafos(t2))));
  }
};

#undef  COMPOSITE_TRAFO




// Traits specialization for Composite:

template<typename TRAFOS>
struct Traits<Composite<TRAFOS> >
{
  typedef Composite<TRAFOS> TRAFO;

  static const int N_RANK=TRAFO::N_RANK;
  
  typedef typename TRAFO::TrafoTypes     TrafoTypes    ;
  typedef typename TRAFO::StateVectorLow StateVectorLow;

  static void transform(const TRAFO& trafo, const StateVectorLow& in, StateVectorLow& out) {trafo.transform(in,out);}

  static const TrafoTypes& getTrafos(const TRAFO& trafo) {return trafo.getTrafos();}

};





} // transformation


} // quantumdata

#endif // _TRANSFORMATION_IMPL_H
