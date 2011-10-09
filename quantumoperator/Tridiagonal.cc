#include "Tridiagonal.h"



namespace quantumoperator {

template<>
const mpl::int_<1> Tridiagonal<1>::_1_ =mpl::int_<1>();

template<>
const Tridiagonal<1>::Diagonal Tridiagonal<1>::empty=Tridiagonal<1>::Diagonal();


namespace {

bool Compatible(const TTD_CARRAY(1)& a, size_t diffA, const TTD_CARRAY(1)& b, size_t diffB=0)
{
  if (!a.size() || !b.size() || a.size()-diffA==b.size()-diffB) return true;
  return false;
}

} 


template<>
Tridiagonal<1>&
Tridiagonal<1>::furnishWithFreqs(const Diagonal& mainDiagonal, mpl::int_<1>)
{
  const size_t k=differences_(0), diagonalSize=getTotalDimension()-k;
  if (diagonalSize>0) {
    freqs_(1).resize(diagonalSize);
    for (int n=0; n<diagonalSize; n++)
      freqs_(1)(n)=mainDiagonal(n+k)-mainDiagonal(n);
    freqs_(2).resize(diagonalSize);
    freqs_(2)=-freqs_(1);
  }
  return *this;
}


template<>
Tridiagonal<1>::Tridiagonal(const Diagonal& zero, size_t k, const Diagonal& minus, const Diagonal& plus, bool toFreqs, mpl::int_<1>)
  : Base(std::max(std::max(size_t(zero.size()),minus.size()+k),plus.size()+k)),
    // funnily enough, Array::size() returns an int ...
    diagonals_(blitzplusplus::TOA_DeepCopy(),
	       toFreqs ? empty : zero,
	       minus,
	       plus),
    differences_(k),
    tCurrent_(0),
    freqs_()
{
  if (toFreqs) furnishWithFreqs(Diagonal(-zero));
  // The minus sign here assures that the dynamics in the two cases of toFreqs are compatible.
  // Consistency check follows:
  if (
      !Compatible(minus,k,zero) || !Compatible(plus,k,zero) || !Compatible(minus,k,plus,k) || 
      ((minus.size() || plus.size()) && !k)
      )
    throw TridiagonalConsistencyErrorException();

}


const Tridiagonal<1> identity(size_t dim)
{
  Tridiagonal<1>::Diagonal diagonal(dim);
  diagonal=1.;
  return Tridiagonal<1>(diagonal);
}




size_t details::binOp1(size_t otherDifference, size_t& difference)
// Declared in Tridiagonal.tcc
{
  if (!difference) return otherDifference;
  else if (!otherDifference || difference==otherDifference) return difference;
  else throw TridiagonalStructureMismatchException();
}



} // quantumoperator


