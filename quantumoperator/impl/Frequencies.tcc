// -*- C++ -*-
#ifndef _FREQUENCIES_IMPL_H
#define _FREQUENCIES_IMPL_H


namespace quantumoperator {


template<int RANK> template<int RANK2>
Frequencies<RANK>::Frequencies(const Frequencies<RANK2>& f1, const Frequencies<RANK-RANK2>& f2)
  : freqs_(blitzplusplus::TOA_ShallowCopy(),details::directDiagonals<RANK2,RANK-RANK2,false>(f1.get(),f2.get()))
{
}


template<int RANK>
Tridiagonal<RANK>& Tridiagonal<RANK>::propagate(const Frequencies<RANK>& f, double t)
{
  const Diagonals& freqs(f.get());
  // Diagonal (i=0) left out
  if (double dt=t-tCurrent_)
    for (int i=1; i<LENGTH; i++)
      if (diagonals_(i).size() && freqs(i).size())
	diagonals_(i)*=exp(dt*freqs(i));
  tCurrent_=t;
  return *this;
}


template<int RANK>
std::ostream& operator<<(std::ostream& os, const Frequencies<RANK>& freqs)
{
  using namespace std;
  os<<freqs.get()<<endl;
  return os;
}


} // quantumoperator


#endif // _FREQUENCIES_IMPL_H
