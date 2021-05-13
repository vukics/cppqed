// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#if !BOOST_PP_IS_ITERATING

#ifndef CPPQEDCORE_QUANTUMOPERATOR_TRIDIAGONAL_TCC_INCLUDED
#define CPPQEDCORE_QUANTUMOPERATOR_TRIDIAGONAL_TCC_INCLUDED

#include "Tridiagonal.h"

#include "ComplexArrayExtensions.h"

#include <boost/lambda/lambda.hpp>
#include <boost/range/algorithm_ext/for_each.hpp>

#include <boost/preprocessor/iteration/iterate.hpp>
#include <boost/preprocessor/repetition/enum.hpp>

#include <boost/mpl/vector_c.hpp>

#include <algorithm>


namespace mpl=boost::mpl;



namespace quantumoperator {


namespace details {


template<bool IS_MULTIPLICATION, int RANK1, int RANK2> // if not multiply then sum
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

      d.reference(doDirect<IS_MULTIPLICATION,RANK1,RANK2>(d1,d2));

    }

  return res;

} // NEEDS_WORK this is now runtime, but (a smaller) part of it could be done compile-time. Ugh. Ouch.


} // details


template<int RANK> template<int RANK2>
Tridiagonal<RANK>::Tridiagonal(const Tridiagonal<RANK2>& t1, const Tridiagonal<RANK-RANK2>& t2)
  : Base(blitzplusplus::concatenateTinies(t1.getDimensions(),t2.getDimensions())),
    diagonals_(blitzplusplus::ShallowCopy(),details::directDiagonals<true,RANK2,RANK-RANK2>(t1.get(),t2.get())),
    differences_(blitzplusplus::concatenateTinies(t1.getDifferences(),t2.getDifferences())),
    tCurrent_(t1.getTime()),
    freqs_(blitzplusplus::ShallowCopy(),details::directDiagonals<false,RANK2,RANK-RANK2>(t1.getFreqs(),t2.getFreqs()))
{
  if (t1.getTime()!=t2.getTime()) throw std::runtime_error("Tridiagonal time mismatch");
}


template<int RANK>
Tridiagonal<RANK>&
Tridiagonal<RANK>::operator*=(dcomp dc)
{
  boost::for_each(diagonals_,boost::lambda::_1*=dc); 
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


template<int RANK>
Tridiagonal<RANK>&
Tridiagonal<RANK>::operator+=(const Tridiagonal& tridiag)
{
  boost::for_each(tridiag.differences_,differences_,[](size_t otherDifference, size_t& difference) {
    if (!difference) difference=otherDifference;
    else if (otherDifference && difference!=otherDifference) throw std::invalid_argument("Tridiagonal structure mismatch is binOp1");
  });
  
  boost::for_each(tridiag.diagonals_,diagonals_,[](const typename Tridiagonal<RANK>::Diagonal& from, typename Tridiagonal<RANK>::Diagonal& to) {
    if (from.size()) {
      if (!to.size()) {
        to.resize(from.shape());
        to=from;
      }
      else to+=from; // This will check for the compatibility of shapes
    }
  });
  
  boost::for_each(tridiag.freqs_,freqs_,[](const typename Tridiagonal<RANK>::Diagonal& from, typename Tridiagonal<RANK>::Diagonal& to) {
    if (from.size()) {
      if (to.size() && all(to!=from)) throw std::runtime_error("Tridiagonal structure mismatch is operator+=");
      else {
        to.resize(from.shape());
        to=from;
      }
    }
  });

  return *this;

}



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
Tridiagonal<RANK>
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
  
  mpl::for_each<mpl::range_c<int,0,RANK>>(FillRangesHelper(res,ubound,differences_));
  return res;
}



template<int RANK>
void Tridiagonal<RANK>::apply(const StateVectorLow& psi, StateVectorLow& dpsidt) const
{
  static const int step=tmptools::Power_v<3,RANK-1>;

  using mpl::vector_c; using tmptools::integral_c;

  const Ranges ranges(fillRanges(psi.ubound()));

  doApply<0*step,vector_c<int,0>,vector_c<int,0>,vector_c<int,0> >(integral_c<RANK-1>(),ranges,psi,dpsidt);
  doApply<1*step,vector_c<int,2>,vector_c<int,1>,vector_c<int,1> >(integral_c<RANK-1>(),ranges,psi,dpsidt);
  doApply<2*step,vector_c<int,1>,vector_c<int,1>,vector_c<int,2> >(integral_c<RANK-1>(),ranges,psi,dpsidt);

}



template<int RANK> template<int START, typename V_DPSIDT, typename V_A, typename V_PSI, int REMAINING>
void Tridiagonal<RANK>::doApply(tmptools::integral_c<REMAINING>,
                                const Ranges& ranges, const StateVectorLow& psi, StateVectorLow& dpsidt) const
{
  static const int step=tmptools::Power_v<3,REMAINING-1>;

  using mpl::push_back, tmptools::integral_c;

  doApply<START+0*step,
          typename push_back<V_DPSIDT,integral_c<0> >::type,
          typename push_back<V_A     ,integral_c<0> >::type,
          typename push_back<V_PSI   ,integral_c<0> >::type>
    (integral_c<REMAINING-1>(),ranges,psi,dpsidt);

  doApply<START+1*step,
          typename push_back<V_DPSIDT,integral_c<2> >::type,
          typename push_back<V_A     ,integral_c<1> >::type,
          typename push_back<V_PSI   ,integral_c<1> >::type>
    (integral_c<REMAINING-1>(),ranges,psi,dpsidt);

  doApply<START+2*step,
          typename push_back<V_DPSIDT,integral_c<1> >::type,
          typename push_back<V_A     ,integral_c<1> >::type,
          typename push_back<V_PSI   ,integral_c<2> >::type>
    (integral_c<REMAINING-1>(),ranges,psi,dpsidt);

}


} // quantumoperator


#define BOOST_PP_ITERATION_LIMITS (1,QUANTUMOPERATOR_TRIDIAGONAL_MAX_RANK)
#define BOOST_PP_FILENAME_1 "Tridiagonal.tcc"

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


#endif // CPPQEDCORE_QUANTUMOPERATOR_TRIDIAGONAL_TCC_INCLUDED


#else  // BOOST_PP_IS_ITERATING


#define rank BOOST_PP_ITERATION()

#define AT_helper(V,i) ranges(i)(at_c<V,i>::type::value)

#define DPSIDT_print(z,m,data) AT_helper(V_DPSIDT,m)
#define      A_print(z,m,data) AT_helper(V_A     ,m)
#define    PSI_print(z,m,data) AT_helper(V_PSI   ,m)


template<> template<int START, typename V_DPSIDT, typename V_A, typename V_PSI>
void quantumoperator::Tridiagonal<rank>::doApply(tmptools::integral_c<0>,
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
