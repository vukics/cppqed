// -*- C++ -*-
#if !BOOST_PP_IS_ITERATING

#ifndef QUANTUMOPERATOR_IMPL_TRIDIAGONAL_TCC_INCLUDED
#define QUANTUMOPERATOR_IMPL_TRIDIAGONAL_TCC_INCLUDED

#include "Tridiagonal.h"

#include "impl/ComplexArrayExtensions.tcc"

#include <boost/lambda/lambda.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/range/algorithm_ext/for_each.hpp>

#include <boost/preprocessor/iteration/iterate.hpp>
#include <boost/preprocessor/repetition/enum.hpp>

#include <algorithm>



namespace quantumoperator {


namespace bll=boost::lambda;



template<int RANK> template<int RANK2>
Tridiagonal<RANK>::Tridiagonal(const Tridiagonal<RANK2>& t1, const Tridiagonal<RANK-RANK2>& t2)
  : Base(blitzplusplus::concatenateTinies(t1.getDimensions(),t2.getDimensions())),
    diagonals_(blitzplusplus::ShallowCopy(),details::directDiagonals<RANK2,RANK-RANK2,true>(t1.get(),t2.get())),
    differences_(blitzplusplus::concatenateTinies(t1.getDifferences(),t2.getDifferences())),
    tCurrent_(t1.getTime()),
    freqs_(blitzplusplus::ShallowCopy(),details::directDiagonals<RANK2,RANK-RANK2,false>(t1.getFreqs(),t2.getFreqs()))
{
  if (t1.getTime()!=t2.getTime()) throw TridiagonalTimeMismatchException();
}


template<int RANK>
Tridiagonal<RANK>&
Tridiagonal<RANK>::operator*=(const dcomp& dc)
{
  boost::for_each(diagonals_,bll::_1*=dc); 
  return *this;
}




template<int RANK>
const Tridiagonal<RANK> Tridiagonal<RANK>::hermitianConjugate() const
{
  struct helper {

    static int transposeIndex(int ind)
    {
      int res=0;
      // std::cout<<ind<<' ';
      for (int k=RANK-1; k>=0; k--) {
	int threeToK=int(floor(pow(3.,k)));
	int ik=int(floor(ind/double(threeToK)));
	res+=((3-ik)%3)*threeToK;
	ind-=ik*threeToK;
      }
      // std::cout<<res<<std::endl;
      return res;
    }

  };

  Tridiagonal res(*this,Diagonals(),differences_,tCurrent_,freqs_);

  for (int ind=0; ind<LENGTH; ind++)
    res.diagonals_(helper::transposeIndex(ind)).reference(Diagonal(conj(diagonals_(ind))));
  // Here a new copy of the conjugated diagonal gets created & referenced.

  // freqs_ at any moment stores the correct transposed values as well

  return res;

}


namespace details {

void binOp1(size_t otherDifference, size_t& difference);

} // details


template<int RANK>
Tridiagonal<RANK>&
Tridiagonal<RANK>::operator+=(const Tridiagonal& tridiag)
{
  struct helper
  {
    static void
    doIt1(const typename Tridiagonal<RANK>::Diagonal& from, typename Tridiagonal<RANK>::Diagonal& to)
    {
      if (from.size()) {
	if (!to.size()) {
	  to.resize(from.shape());
	  to=from;
	}
	else to+=from; // This will check for the compatibility of shapes
      }
    }

    static void
    doIt2(const typename Tridiagonal<RANK>::Diagonal& from, typename Tridiagonal<RANK>::Diagonal& to)
    {
      if (from.size()) {
	if (to.size() && all(to!=from)) throw TridiagonalStructureMismatchException();
	else {
	  to.resize(from.shape());
	  to=from;
	}
      }
    }

  };

  boost::for_each(tridiag.differences_,differences_,details::binOp1);
  boost::for_each(tridiag.  diagonals_,  diagonals_, helper::doIt1 );
  boost::for_each(tridiag.      freqs_,      freqs_, helper::doIt2 );

  return *this;

}



namespace details {



template<int RANK1, int RANK2, bool MULT> // if not multiply then sum
const typename Tridiagonal<RANK1+RANK2>::Diagonals
directDiagonals(const typename Tridiagonal<RANK1>::Diagonals& ds1,
		const typename Tridiagonal<RANK2>::Diagonals& ds2) 
{
  using namespace blitzplusplus; // implies blitz
  using namespace linalg       ;

  typedef Tridiagonal<RANK1+RANK2> TridiagonalRes;
  typedef Tridiagonal<RANK1> Tridiagonal1;
  typedef Tridiagonal<RANK2> Tridiagonal2;


  typename TridiagonalRes::Diagonals res;
  size_t length=ds2.numElements();

  for (size_t i=0; i<ds1.numElements(); ++i) for (size_t j=0; j<length; ++j) {
      const typename Tridiagonal1::Diagonal& d1=ds1(i);
      const typename Tridiagonal2::Diagonal& d2=ds2(j);
      typename TridiagonalRes::Diagonal& d=res(i*length+j);

      d.reference(doDirect(d1,d2,mpl::bool_<MULT>()));

    }

  return res;

} // NEEDS_WORK this is now runtime, but (a smaller) part of it could be done compile-time. Ugh. Ouch.


} // details



template<int RANK>
Tridiagonal<RANK>& Tridiagonal<RANK>::propagate(double t)
{
  // Diagonal (i=0) left out
  if (double dt=t-tCurrent_)
    for (int i=1; i<LENGTH; i++)
      if (diagonals_(i).size() && freqs_(i).size())
	diagonals_(i)*=exp(dt*freqs_(i));
  tCurrent_=t;
  return *this;
}


template<int RANK>
const Tridiagonal<RANK>
furnishWithFreqs(const Tridiagonal<RANK>& tridiag, const typename Tridiagonal<RANK>::Diagonal& mainDiagonal)
{
  Tridiagonal<RANK> res(tridiag);
  return res.furnishWithFreqs(mainDiagonal);
}


template<int RANK>
std::ostream& operator<<(std::ostream& os, const Tridiagonal<RANK>& tridiag)
{
  using namespace std;
  os<<tridiag.get()<<endl<<tridiag.getFreqs()<<endl<<tridiag.getDifferences()<<endl<<tridiag.getDimensions()<<endl;
  return os;
}


} // quantumoperator


#ifndef DO_CONSIDER_EXPLICITLY_SPECIALIZED_TRIDIAGONAL_APPLIES


namespace quantumoperator {


template<int RANK> template<typename ICW>
void Tridiagonal<RANK>::FillRangesHelper::operator()(ICW)
{
  using blitz::Range;

  static const int i=ICW::value;
  
  ranges_(i)(0)=Range(0,ubound_(i));
  ranges_(i)(2)=((
		  ranges_(i)(1)=Range(0,ubound_(i)-int(k_(i)))
		  )+int(k_(i)));

}



template<int RANK>
const typename Tridiagonal<RANK>::Ranges Tridiagonal<RANK>::fillRanges(const typename StateVectorLow::T_index& ubound) const
{
  Ranges res;
  
  mpl::for_each<tmptools::Ordinals<RANK> >(FillRangesHelper(res,ubound,differences_));
  return res;
}



template<int RANK>
void Tridiagonal<RANK>::apply(const StateVectorLow& psi, StateVectorLow& dpsidt) const
{
  static const int step=tmptools::Power<3,RANK-1>::value;

  using tmptools::Vector; using mpl::int_;

  const Ranges ranges(fillRanges(psi.ubound()));

  doApply<0*step,Vector<0>,Vector<0>,Vector<0> >(int_<RANK-1>(),ranges,psi,dpsidt);
  doApply<1*step,Vector<2>,Vector<1>,Vector<1> >(int_<RANK-1>(),ranges,psi,dpsidt);
  doApply<2*step,Vector<1>,Vector<1>,Vector<2> >(int_<RANK-1>(),ranges,psi,dpsidt);

}



template<int RANK> template<int START, typename V_DPSIDT, typename V_A, typename V_PSI, int REMAINING>
void Tridiagonal<RANK>::doApply(mpl::int_<REMAINING>,
				const Ranges& ranges, const StateVectorLow& psi, StateVectorLow& dpsidt) const
{
  static const int step=tmptools::Power<3,REMAINING-1>::value;

  using mpl::push_back; using mpl::int_;

  doApply<START+0*step,
	  typename push_back<V_DPSIDT,int_<0> >::type,
	  typename push_back<V_A     ,int_<0> >::type,
	  typename push_back<V_PSI   ,int_<0> >::type>
    (int_<REMAINING-1>(),ranges,psi,dpsidt);

  doApply<START+1*step,
	  typename push_back<V_DPSIDT,int_<2> >::type,
	  typename push_back<V_A     ,int_<1> >::type,
	  typename push_back<V_PSI   ,int_<1> >::type>
    (int_<REMAINING-1>(),ranges,psi,dpsidt);

  doApply<START+2*step,
	  typename push_back<V_DPSIDT,int_<1> >::type,
	  typename push_back<V_A     ,int_<1> >::type,
	  typename push_back<V_PSI   ,int_<2> >::type>
    (int_<REMAINING-1>(),ranges,psi,dpsidt);

}


} // quantumoperator


#define BOOST_PP_ITERATION_LIMITS (1,QUANTUMOPERATOR_TRIDIAGONAL_MAX_RANK)
#define BOOST_PP_FILENAME_1 "impl/Tridiagonal.tcc"

#include BOOST_PP_ITERATE()

#undef BOOST_PP_FILENAME_1
#undef BOOST_PP_ITERATION_LIMITS


#else // DO_CONSIDER_EXPLICITLY_SPECIALIZED_TRIDIAGONAL_APPLIES

namespace quantumoperator {


template<>
void Tridiagonal<1>::apply(const StateVectorLow&, StateVectorLow&) const;

template<>
void Tridiagonal<2>::apply(const StateVectorLow&, StateVectorLow&) const;

template<>
void Tridiagonal<3>::apply(const StateVectorLow&, StateVectorLow&) const;

template<>
void Tridiagonal<4>::apply(const StateVectorLow&, StateVectorLow&) const;

template<>
void Tridiagonal<5>::apply(const StateVectorLow&, StateVectorLow&) const;

template<>
void Tridiagonal<6>::apply(const StateVectorLow&, StateVectorLow&) const;

template<>
void Tridiagonal<7>::apply(const StateVectorLow&, StateVectorLow&) const;

template<>
void Tridiagonal<8>::apply(const StateVectorLow&, StateVectorLow&) const;

template<>
void Tridiagonal<9>::apply(const StateVectorLow&, StateVectorLow&) const;


} // quantumoperator


#endif // DO_CONSIDER_EXPLICITLY_SPECIALIZED_TRIDIAGONAL_APPLIES


#endif // QUANTUMOPERATOR_IMPL_TRIDIAGONAL_TCC_INCLUDED


#else  // BOOST_PP_IS_ITERATING


#define rank BOOST_PP_ITERATION()

#define AT_helper(V,i) ranges(i)(at_c<V,i>::type::value)

#define DPSIDT_print(z,m,data) AT_helper(V_DPSIDT,m)
#define      A_print(z,m,data) AT_helper(V_A     ,m)
#define    PSI_print(z,m,data) AT_helper(V_PSI   ,m)


template<> template<int START, typename V_DPSIDT, typename V_A, typename V_PSI>
void quantumoperator::Tridiagonal<rank>::doApply(mpl::int_<0>,
                                                 const Ranges& ranges, const StateVectorLow& psi, StateVectorLow& dpsidt) const
{
  using mpl::at_c;

  if (diagonals_(START).size()) 
    dpsidt             (BOOST_PP_ENUM(rank,DPSIDT_print,~))
      +=
      diagonals_(START)(BOOST_PP_ENUM(rank,     A_print,~))
      *
      psi              (BOOST_PP_ENUM(rank,   PSI_print,~))
      ;

}

#undef    PSI_print
#undef      A_print
#undef DPSIDT_print

#undef AT_helper

#undef rank


#endif // BOOST_PP_IS_ITERATING
