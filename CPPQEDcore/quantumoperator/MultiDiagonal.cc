// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "MultiDiagonal.h"


auto quantumoperator::operator|(const quantumoperator::MultiDiagonal<1>& a, const quantumoperator::MultiDiagonal<1>& b) -> MultiDiagonal<1>
{
  using Offsets=MultiDiagonal<1>::Offsets;
  using Diagonal=MultiDiagonal<1>::Diagonal;

  size_t dimension=calculateAndCheckDimensions(a)[0];
  if (dimension!=calculateAndCheckDimensions(b)[0]) throw std::runtime_error("Mismatch in MultiDiagonal::compose dimensions");

  MultiDiagonal<1> res;

  for (const auto& [aidx,adti] : a.diagonals) for (const auto& [bidx,bdti] : b.diagonals) for (const auto& [aoffset,adiag] : adti) for (const auto& [boffset,bdiag] : bdti) {

    auto diagInit = [&] (ptrdiff_t i, ptrdiff_t j, size_t start=0) {
      return [a=adiag.dataView, b=bdiag.dataView, i=i, j=j, start=start] (size_t e) {
        auto res{zeroInit<1>(e)}; for (size_t k=start; k<e; ++k) res[k]=a[k+i]*b[k+j]; return res;
      };
    };

    auto emplaceDiagonal = [&] (size_t idx, size_t offset, size_t extent, ptrdiff_t i, ptrdiff_t j, size_t start=0) {
      res.diagonals[idx].emplace(Offsets{offset}, Diagonal{{extent}, diagInit(i,j,start)});
    };

    // both upper
    if (size_t n=aoffset[0], m=boffset[0]; aidx==1 && bidx==1) {
      if (m+n<dimension) emplaceDiagonal( 1, n+m, adiag.extents[0]-m, 0, n ) ;
    }
    // both lower
    else if (aidx==0 && bidx==0) {
      if (m+n<dimension) emplaceDiagonal( 0, n+m, adiag.extents[0]-m, m, 0 ) ;
    }
    // upper with lower
    else if (aidx==1 && bidx==0) {
      if (n<m) // the result is lower
        emplaceDiagonal( 0, m-n, bdiag.extents[0], m-n, 0 ) ;
      else // the result is upper
        emplaceDiagonal( 1, n-m, adiag.extents[0], 0, n-m ) ;
    }
    // lower with upper
    else {
      if (m<n) // the result is lower
        emplaceDiagonal( 0, n-m, adiag.extents[0]+m, -m, -m, m ) ;
      else // the result is upper
        emplaceDiagonal( 1, m-n, bdiag.extents[0]+n, -n, -n, n ) ;
    }
  }

  return res;
}



/*
namespace quantumoperator {

template<>
const Tridiagonal<1>::Diagonal Tridiagonal<1>::empty=Tridiagonal<1>::Diagonal();


namespace {

bool Compatible(const CArray<1>& a, size_t diffA, const CArray<1>& b, size_t diffB=0)
{
  if (!a.size() || !b.size() || a.size()-diffA==b.size()-diffB) return true;
  return false;
}

} 


template<> template<>
Tridiagonal<1>&
Tridiagonal<1>::furnishWithFreqs(const Diagonal& mainDiagonal)
{
  const size_t k=differences_(0), diagonalSize=getTotalDimension()-k;
  if (diagonalSize>0) {
    freqs_(1).resize(diagonalSize);
    for (int n=0; n<diagonalSize; n++)
      freqs_(1)(n)=mainDiagonal(n+k)-mainDiagonal(n);
    freqs_(2).reference(Diagonal(-freqs_(1)));
  }
  return *this;
}


template<> template<>
Tridiagonal<1>::Tridiagonal(const Diagonal& zero, size_t k, const Diagonal& minus, const Diagonal& plus, bool toFreqs)
  : Base(std::max(std::max(size_t(zero.size()),minus.size()+k),plus.size()+k)),
    // funnily enough, Array::size() returns an int ...
    diagonals_(blitzplusplus::DeepCopy(),
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
      !Compatible(zero,k,minus) || !Compatible(zero,k,plus) || !Compatible(minus,0,plus,0) || 
      ((minus.size() || plus.size()) && !k)
      )
    throw std::invalid_argument("Tridiagonal consistency error");

}


const Tridiagonal<1> identity(size_t dim)
{
  Tridiagonal<1>::Diagonal diagonal(dim);
  diagonal=1.;
  return Tridiagonal<1>(diagonal);
}



} // quantumoperator


*/
