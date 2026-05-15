// Copyright András Vukics 2006–2026. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "MultiDiagonal.h"


/// @brief Rank-1 identity operator: a single main diagonal (`Index=1`, `Offsets={0}`) filled with 1.
quantumoperator::MultiDiagonal<1> quantumoperator::multidiagonal::identity(size_t dim)
{
  MultiDiagonal<1> res;
  res.diagonals[1].emplace(
    MultiDiagonal<1>::Offsets{0},
    MultiDiagonal<1>::Diagonal{ multiarray::StorageType<dcomp>(dim,1.) } );
  return res;
}


/// @brief Rank-1 composition \f$(A|B)_{nk} = \sum_m A_{nm} B_{mk}\f$.
///
/// Iterates over all pairs of diagonals from @p a and @p b. For each pair the result offset
/// and direction (upper/lower) are determined by four cases based on the `Index` bits:
///
/// | a         | b         | condition | result    | result offset |
/// |-----------|-----------|-----------|-----------|---------------|
/// | upper (1) | upper (1) | n+m < dim | upper     | n+m           |
/// | lower (0) | lower (0) | n+m < dim | lower     | n+m           |
/// | upper (1) | lower (0) | n < m     | lower     | m−n           |
/// | upper (1) | lower (0) | n ≥ m     | upper     | n−m           |
/// | lower (0) | upper (1) | m < n     | lower     | n−m           |
/// | lower (0) | upper (1) | m ≥ n     | upper     | m−n           |
///
/// The result diagonal is built by pointwise multiplication of the two source dataViews
/// on the valid index range, via an initializer lambda passed to `MultiArray`.
///
/// @note `operator*` has higher precedence than `operator|` in C++; use parentheses in mixed expressions.
/// @note Higher-rank composition requires generalization of the four upper/lower cases above.
auto quantumoperator::operator|(const quantumoperator::MultiDiagonal<1>& a, const quantumoperator::MultiDiagonal<1>& b) -> MultiDiagonal<1>
{
  using Offsets=MultiDiagonal<1>::Offsets;
  using Diagonal=MultiDiagonal<1>::Diagonal;

  size_t dimension=calculateAndCheckDimensions(a)[0];
  if (dimension!=calculateAndCheckDimensions(b)[0]) throw std::runtime_error("Mismatch in MultiDiagonal::compose dimensions");

  MultiDiagonal<1> res;

  auto diagInit = [&] (ptrdiff_t i, ptrdiff_t j, size_t start=0) {
    return [a=a.diagonals.begin()->second.begin()->second.dataView,
            b=b.diagonals.begin()->second.begin()->second.dataView,
            i=i, j=j, start=start] (size_t e) {
      auto res{multiarray::zeroInit<dcomp>(e)}; for (size_t k=start; k<e; ++k) res[k]=a[k+i]*b[k+j]; return res;
    };
  };

  auto emplaceDiagonal = [&] (size_t idx, size_t offset, size_t extent, ptrdiff_t i, ptrdiff_t j, size_t start=0) {
    res.diagonals[idx].emplace(Offsets{offset}, Diagonal{{extent}, diagInit(i,j,start)});
  };

  for (const auto& [aidx,adti] : a.diagonals) for (const auto& [bidx,bdti] : b.diagonals)
    for (const auto& [aoffset,adiag] : adti) for (const auto& [boffset,bdiag] : bdti) {

    size_t n=aoffset[0], m=boffset[0];

    // both upper
    if (aidx==1 && bidx==1) {
      if (m+n<dimension) emplaceDiagonal( 1, n+m, adiag.extents[0]-m, 0, n );
    }
    // both lower
    else if (aidx==0 && bidx==0) {
      if (m+n<dimension) emplaceDiagonal( 0, n+m, adiag.extents[0]-m, m, 0 );
    }
    // upper with lower
    else if (aidx==1 && bidx==0) {
      if (n<m) // result is lower
        emplaceDiagonal( 0, m-n, bdiag.extents[0], m-n, 0 );
      else     // result is upper
        emplaceDiagonal( 1, n-m, adiag.extents[0], 0, n-m );
    }
    // lower with upper
    else {
      if (m<n) // result is lower
        emplaceDiagonal( 0, n-m, adiag.extents[0]+m, -m, -m, m );
      else     // result is upper
        emplaceDiagonal( 1, m-n, bdiag.extents[0]+n, -n, -n, n );
    }
  }

  return res;
}