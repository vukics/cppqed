// -*- C++ -*-
#ifndef _FREQUENCIES_H
#define _FREQUENCIES_H

// An interesting idea would be to fuse Frequencies into a
// TridiagonalInteractionPicture class which could be derived from
// Tridiagonal. The problem which has prevented me from doing this is
// that I don't know how to lucidly define addition for the embedded
// Frequencies when adding two such classes.

#include "Tridiagonal.h"


namespace quantumoperator {


// const Frequencies<1> zeroFreq();


template<int RANK>
class Frequencies
{  
public:
  typedef typename Tridiagonal<RANK>::Diagonals Diagonals;
  typedef typename Tridiagonal<RANK>::Diagonal  Diagonal ;

  Frequencies(size_t size, mpl::int_<RANK> =mpl::int_<1>());

  Frequencies(size_t K, const Diagonal& minus, mpl::int_<RANK> =mpl::int_<1>());

  Frequencies(const Frequencies& freqs) 
    : freqs_(blitzplusplus::TOA_DeepCopy(),freqs.freqs_) {}

  template<int RANK2>
  Frequencies(const Frequencies<RANK2>&, const Frequencies<RANK-RANK2>&); // Direct product (sum)


  const Diagonals& get() const {return freqs_;}

private:
  Diagonals freqs_;

};


template<int RANK1, int RANK2>
inline
const Frequencies<RANK1+RANK2>
operator*(const Frequencies<RANK1>& t1, const Frequencies<RANK2>& t2)
{
  return Frequencies<RANK1+RANK2>(t1,t2);
}


template<int RANK>
std::ostream& operator<<(std::ostream&, const Frequencies<RANK>&);


} // quantumoperator


#include "impl/Frequencies.tcc"


#endif // _FREQUENCIES_H
