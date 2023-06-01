// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "MultiDiagonal.h"



auto quantumoperator::compose(const quantumoperator::MultiDiagonal<1>& a, const quantumoperator::MultiDiagonal<1> b) -> MultiDiagonal<1>
{
  if (a.dimensions!=b.dimensions) throw std::runtime_error("Mismatch in MultiDiagonal::compose dimensions");

  MultiDiagonal<1> res{a.dimensions};

  for (auto&& [ao,ad] : a.diagonals) for (auto&& [bo,bd] : b.diagonals) {
    int n=ao[0], s=ao[0]+bo[0];
    MultiDiagonal<1>::Offsets composedOffsets{s};
    MultiDiagonal<1>::Diagonal composedDiagonal{ {a.dimensions[0]-s}, ::quantumdata::noInit<1> };
    auto composedIndexer{MultiDiagonal<1>::diagonalIndexer(composedOffsets,composedDiagonal)};
    auto aIndexer{MultiDiagonal<1>::diagonalIndexer(ao,ad)};
    auto bIndexer{MultiDiagonal<1>::diagonalIndexer(bo,bd)};
    for (size_t k=0; k<a.dimensions[0]-s; ++k) composedIndexer({k})=aIndexer({k})*bIndexer({k+n});
    try { for (auto&& [de,cde] : std::views::zip(res.diagonals.at(composedOffsets).mutableView().dataView,composedDiagonal.dataView)) de+=cde; }
    catch (const std::out_of_range&) { res.diagonals.insert({composedOffsets,std::move(composedDiagonal)}); }
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
