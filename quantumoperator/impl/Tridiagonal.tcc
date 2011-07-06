// -*- C++ -*-
#ifndef _TRIDIAGONAL_IMPL_H
#define _TRIDIAGONAL_IMPL_H

#include "ComplexArrayExtensions.h"

#include <algorithm>



namespace quantumoperator {


namespace bll=boost::lambda;



template<int RANK> template<int RANK2>
Tridiagonal<RANK>::Tridiagonal(const Tridiagonal<RANK2>& t1, const Tridiagonal<RANK-RANK2>& t2)
  : Base(blitzplusplus::concatenateTinies(t1.getDimensions(),t2.getDimensions())),
    diagonals_(blitzplusplus::TOA_ShallowCopy(),details::directDiagonals<RANK2,RANK-RANK2,true>(t1.get(),t2.get())),
    differences_(blitzplusplus::concatenateTinies(t1.getDifferences(),t2.getDifferences())),
    tCurrent_(0),
    freqs_(blitzplusplus::TOA_ShallowCopy(),details::directDiagonals<RANK2,RANK-RANK2,false>(t1.getFreqs(),t2.getFreqs()))
{
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

  // freqs_ at any moment stores the correct transposed values as well

  return res;

}


namespace details {

size_t binOp1(size_t otherDifference, size_t& difference);

} // details


template<int RANK>
Tridiagonal<RANK>&
Tridiagonal<RANK>::operator+=(const Tridiagonal& tridiag)
{
  struct helper
  {
    static typename Tridiagonal<RANK>::Diagonal& 
    doIt1(const typename Tridiagonal<RANK>::Diagonal& from, typename Tridiagonal<RANK>::Diagonal& to)
    {
      if (from.size()) {
	if (!to.size()) {
	  to.resize(from.shape());
	  to=from;
	}
	else to+=from; // This will check for the compatibility of shapes
      }
      return to;
    }

    static typename Tridiagonal<RANK>::Diagonal& 
    doIt2(const typename Tridiagonal<RANK>::Diagonal& from, typename Tridiagonal<RANK>::Diagonal& to)
    {
      if (from.size()) {
	if (to.size() && all(to!=from)) throw TridiagonalStructureMismatchException();
	else {
	  to.resize(from.shape());
	  to=from;
	}
      }
      return to;
    }

  };

  boost::transform(tridiag.differences_,differences_.begin(),differences_.begin(),details::binOp1);
  boost::transform(tridiag.  diagonals_,  diagonals_.begin(),  diagonals_.begin(), helper::doIt1 );
  boost::transform(tridiag.      freqs_,      freqs_.begin(),      freqs_.begin(), helper::doIt2 );

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
  size_t length=ds2.numElements;

  for (size_t i=0; i<ds1.numElements; ++i) for (size_t j=0; j<length; ++j) {
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



#endif // _TRIDIAGONAL_IMPL_H
