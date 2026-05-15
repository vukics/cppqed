// Copyright András Vukics 2022–2026. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "MultiArray.h"

#include <Eigen/Dense>

#include <boost/serialization/complex.hpp>

#include <complex>

/// @brief Double-precision complex number. Named to resemble a built-in type rather than a class.
typedef std::complex<double> dcomp;

/// @brief Brings `operator""i` into scope for imaginary literals (e.g. `1.0i`, `0.5i`).
using std::complex_literals::operator""i;

/// @brief Returns `true` if either the real or imaginary part of @p c is non-zero.
inline bool isNonZero(dcomp c) {return bool(real(c)) || bool(imag(c));}

/// @brief Returns `true` if the real part of @p c is non-zero.
inline bool hasRealPart(dcomp c) {return bool(real(c));}
/// @brief Returns `true` if the imaginary part of @p c is non-zero.
inline bool hasImagPart(dcomp c) {return bool(imag(c));}

/// @brief Comparator by absolute value, for use with sorting algorithms.
inline bool  absCompare(dcomp c1, dcomp c2) {return  abs(c1)< abs(c2);}
/// @brief Comparator by real part, for use with sorting algorithms.
inline bool realCompare(dcomp c1, dcomp c2) {return real(c1)<real(c2);}

/// @brief Generic square: `x*x`. Return type is deduced; `noexcept` and SFINAE are propagated.
/// @see https://stackoverflow.com/a/64849248/1171157
auto sqr(auto&& x) // return type is non-reference or trailing
noexcept(noexcept(x*x)) // propagate noexcept
-> decltype(x*x) // enable SFINAE
{ return x * x; }

/// @brief Relative deviation between two values: `|a-b| / (|a|+|b|)`. Returns 0 if both are zero.
double relativeDeviation(const auto& a, const auto& b) {return std::abs(a-b)/(std::abs(a)+std::abs(b));}

/// @brief Squared absolute value of a complex number: `|v|^2 = re^2 + im^2`.
inline double sqrAbs(dcomp v) {return sqr(std::abs(v));}


namespace boost::json {
/// @brief Boost.JSON serialization for `std::complex<T>`: emits a two-element JSON array `[re, im]`.
    template<class T>
    void tag_invoke(value_from_tag, value& jv, const std::complex<T>& c) {
        jv = {c.real(), c.imag()};
    }

/// @brief Boost.JSON deserialization for `std::complex<T>`: reads from a two-element JSON array `[re, im]`.
    template<class T>
    std::complex<T> tag_invoke(value_to_tag<std::complex<T>>, const value& jv) {
        return {jv.at(0), jv.at(1)};
    }
} // boost::json

/// @brief Dense complex matrix type (Eigen dynamic-size).
using CMatrix=Eigen::MatrixX<dcomp>;
/// @brief Dense complex vector type (Eigen dynamic-size).
using CVector=Eigen::VectorX<dcomp>;


namespace cppqedutils {

/// @brief Concept satisfied by scalar numeric types usable as operator coefficients: any type convertible to `double` or `dcomp`.
/// @todo Generalize to arbitrary floating-point or arbitrary-precision types and their complex counterparts.
template <typename N> concept scalar = std::convertible_to<N,double> || std::convertible_to<N,dcomp>;

/// @brief Maps a mutable `MultiArray<dcomp,RANK>` to an Eigen column vector (no data copy).
template <size_t RANK>
auto vectorize(MultiArray<dcomp,RANK>& ma)
{
  return Eigen::Map<CVector>{ma.dataStorage().data(),ma.dataView.size()};
}

/// @brief Maps a const `MultiArray<dcomp,RANK>` to a const Eigen column vector (no data copy).
template <size_t RANK>
auto vectorize(const MultiArray<dcomp,RANK>& ma)
{
  return Eigen::Map<const CVector>{ma.dataView.data(),ma.dataView.size()};
}


/// @brief Extracts the first half of a rank-`2R` extent array, verifying that both halves are equal.
///
/// A `MultiArray<dcomp, 2R>` with equal first and second half extents represents a square matrix
/// (or tensor product of square matrices) of total dimension `prod(extents[0..R-1])`.
/// Throws `std::invalid_argument` in debug mode if the two halves differ.
template <size_t TWO_TIMES_RANK> requires ( !(TWO_TIMES_RANK%2) )
Extents<TWO_TIMES_RANK/2> halveExtents(Extents<TWO_TIMES_RANK> extents)
{
  static constexpr size_t RANK=TWO_TIMES_RANK/2;
  Extents<RANK> res;
  std::ranges::for_each(std::views::iota(0ul,RANK), [&] (size_t i) {
    res[i]=extents[i];
#ifndef   NDEBUG
    if (size_t b=extents[i+RANK]; res[i]!=b)
      throw std::invalid_argument("MultiArray does not have the shape of a matrix: index no. "+std::to_string(i)+": "+std::to_string(res[i])+" vs. "+std::to_string(b));
#endif // NDEBUG
  });
  return res;
}

/// @brief Maps a mutable rank-`2R` `MultiArray<dcomp>` to an Eigen square matrix (no data copy).
/// The matrix dimension is `prod(extents[0..R-1])`; requires equal first and second half extents.
template <size_t TWO_TIMES_RANK>
auto matricize(MultiArray<dcomp,TWO_TIMES_RANK>& ma)
{
  const size_t matrixDim = multiarray::calculateExtent(halveExtents(ma.extents));
  return Eigen::Map<CMatrix>{ma.dataStorage().data(),matrixDim,matrixDim};
}

/// @brief Maps a const rank-`2R` `MultiArray<dcomp>` to a const Eigen square matrix (no data copy).
template <size_t TWO_TIMES_RANK>
auto matricize(const MultiArray<dcomp,TWO_TIMES_RANK>& ma)
{
  const size_t matrixDim = multiarray::calculateExtent(halveExtents(ma.extents));
  return Eigen::Map<const CMatrix>{ma.dataView.data(),matrixDim,matrixDim};
}


/// @brief Tensor (outer) product of two complex multi-arrays.
///
/// Returns a `MultiArray<dcomp, RANK1+RANK2>` with extents `concatenate(m1.extents, m2.extents)`.
/// Element `(i,j)` of the result is `func(m1[i], m2[j])`; defaults to pointwise multiplication.
/// Storage layout: `result[i + stride*j]` where `stride = m1.dataView.size()`.
template <size_t RANK1, size_t RANK2>
auto directProduct(const MultiArray<dcomp,RANK1>& m1, const MultiArray<dcomp,RANK2>& m2,
                   std::function<dcomp(dcomp,dcomp)> func=std::multiplies<dcomp>{} )
{
  MultiArray<dcomp,RANK1+RANK2> res{concatenate(m1.extents,m2.extents)};
  for (size_t stride=m1.dataView.size(), i=0; i<stride; ++i) for (size_t j=0; j<m2.dataView.size(); ++j) res.dataStorage()[i+stride*j]=func(m1.dataView[i],m2.dataView[j]);
  return res;
}


/// @brief In-place complex conjugation of all elements of @p ma.
template <size_t RANK>
MultiArray<dcomp,RANK>& conj(MultiArray<dcomp,RANK>& ma)
{
  for (dcomp& v : ma.dataStorage()) v=conj(v);
  return ma;
}


/// @brief Squared Frobenius norm: \f$\sum_{i}|a_i|^2\f$. Returns 0 for an empty array.
template <size_t RANK>
double frobeniusNormSqr(const MultiArray<dcomp,RANK>& ma) { return std::ranges::fold_left_first( ma.dataView | std::views::transform(sqrAbs), std::plus{} ).value_or(0.); }

/// @brief Frobenius norm: \f$\sqrt{\sum_{i}|a_i|^2}\f$.
template <size_t RANK>
double frobeniusNorm(const MultiArray<dcomp,RANK>& ma) { return sqrt( frobeniusNormSqr(ma) ); }

/// @brief In-place Hermitian conjugation of a rank-`2R` array interpreted as a square matrix.
///
/// Transposes and conjugates simultaneously: swaps `(i,j)` and `(j,i)` while conjugating both,
/// then conjugates the diagonal in place. Operates directly on flat storage.
template <size_t TWO_TIMES_RANK>
void hermitianConjugateSelf(MultiArray<dcomp,TWO_TIMES_RANK>& ma)
{
  const size_t matrixDim = multiarray::calculateExtent(halveExtents(ma.extents));
  for (size_t i=0; i<matrixDim; ++i) for (size_t j=i+1; j<matrixDim; ++j) {
    dcomp
      &u=ma.dataStorage()[i+matrixDim*j],
      &l=ma.dataStorage()[j+matrixDim*i];
    u=conj(u); l=conj(l);
    std::swap(u, l);
  }
}


/// @brief In-place replacement of a rank-`2R` matrix view @p mav with `M + M†` (twice its Hermitian part).
///
/// For each pair `(i,j)` with `j >= i`, sets `mav(i,j) = mav(i,j) + conj(mav(j,i))` and
/// `mav(j,i) = conj(mav(i,j))`. The diagonal becomes `2*re(mav(i,i))`.
/// Used in the Lindblad master equation to enforce Hermiticity after each ODE step.
template <size_t TWO_TIMES_RANK>
void twoTimesRealPartOfSelf(MultiArrayView<dcomp,TWO_TIMES_RANK> mav)
{
  const size_t matrixDim = multiarray::calculateExtent(halveExtents(mav.extents));
  auto _=[&] (size_t i, size_t j) -> dcomp& {return mav.dataView[i+matrixDim*j];};
  for (size_t i=0; i<matrixDim; ++i) {
    _(i,i)=2.*real(_(i,i));
    for (size_t j=i; j<matrixDim; ++j) _(j,i)=conj(_(i,j)=_(i,j)+conj(_(j,i)));
  }
};


} // cppqedutils