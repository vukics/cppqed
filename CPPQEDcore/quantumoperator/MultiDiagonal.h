// Copyright András Vukics 2006–2026. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "MultiArrayComplex.h"

#include <bitset>
#include <map>


/// @brief Comprises modules representing operators of special structure (multidiagonal, sparse) over Hilbert spaces of arbitrary arity.
namespace quantumoperator {

using namespace cppqedutils;

namespace multidiagonal {


/// @brief Returns a flat range of `(index, offsets, diagonal)` triples over all diagonals of @p md.
/// Accepts either a pointer or a reference to a `MultiDiagonal`.
auto range(auto&& md)
{
  if constexpr (std::is_pointer_v<std::decay_t<decltype(md)>>) return range(*md);
  else return md.diagonals | std::views::transform([&](auto&& outer_entry) {
    return outer_entry.second | std::views::transform([&](auto&& inner_entry) {
      return std::tie(outer_entry.first, inner_entry.first, inner_entry.second);
    });
  }) | std::views::join;
}


} // multidiagonal


/// @brief A sparse matrix storing an arbitrary set of diagonals, closed under composition, direct product, and Hermitian conjugation.
///
/// `MultiDiagonal<RANK>` supersedes `Tridiagonal<RANK>` from v2, which was limited to at most three
/// diagonals per axis and was **not closed under composition** — the product of two tridiagonal operators
/// is in general pentadiagonal. `MultiDiagonal` lifts both restrictions:
/// - Any number of diagonals at arbitrary offsets per axis.
/// - `operator|` (composition) produces new diagonals at all pairwise offset sums — the algebra is closed.
/// - `operator*` (direct product) concatenates offset arrays, reflecting tensor product structure exactly.
///
/// This fulfills the design goal stated in the `Tridiagonal` documentation:
/// \f[
///   H/i = \bigotimes_{m=0}^{M-1}
///     \sum_{i_m \in \mathbb{K}_m} \sum_{n_m=0}^{N_m-1-i_m}
///     \alpha^{i_m}_{m,n}\, |n_m+i_m\rangle\langle n_m|
/// \f]
/// with \f$\mathbb{K}_m\f$ an **arbitrary** set of offsets per axis.
///
/// ### Time dependence
/// `MultiDiagonal` is time-independent. Interaction-picture time dependence — diagonal elements
/// of the form \f$\alpha_{m,n} e^{\delta_{m,n} t}\f$ with complex \f$\delta\f$ — is handled by
/// `InteractionPictureDiagonal`, which derives from `MultiDiagonal` and adds a parallel frequency
/// `MultiDiagonal` and `propagate(t)` bookkeeping. This separation keeps `MultiDiagonal` a clean,
/// unconditionally time-independent algebraic type.
///
/// ### Operator conventions
/// - `operator()` applies \f$H/i\f$ (not \f$H\f$) to the state vector, consistent with the framework convention.
/// - \f$H\f$ need not be Hermitian: in QJMC and Master equation contexts it is the full effective
///   non-Hermitian Hamiltonian \f$H_\text{eff} = H_\text{phys} - \frac{i}{2}\sum_k J_k^\dagger J_k\f$.
///
/// ### Design note
/// The `Index` (bitset) + `Offsets` (size_t array) two-level map is a known improvement target:
/// replacing with a single `std::map<std::array<ptrdiff_t,RANK>, Diagonal>` of signed offsets
/// would simplify composition arithmetic and eliminate the two-level structure. Deferred to a future refactor.
template <size_t RANK>
struct MultiDiagonal
{
  using Dimensions = Extents<RANK>;
  /// @brief Direction pattern: `true` = upper (or main) diagonal along axis `i`, `false` = lower.
  /// @note `std::bitset` has hash but no comparison; `std::array` has lexicographic comparison by default.
  using Index = std::bitset<RANK>;
  using Offsets = Dimensions;
  using Diagonal = MultiArray<dcomp,RANK>;
  using DiagToIdx = std::map<Offsets,Diagonal>;

  /// @brief Two-level map: direction pattern → offset magnitude → diagonal data.
  using Diagonals = std::unordered_map<Index,DiagToIdx>;

  MultiDiagonal(const MultiDiagonal&) = delete; MultiDiagonal& operator=(const MultiDiagonal&) = delete;
  MultiDiagonal(MultiDiagonal&&) = default; MultiDiagonal& operator=(MultiDiagonal&&) = default;

  explicit MultiDiagonal(auto&&... args) : diagonals{FWD(args)...} {}

  /// @brief Named deep copy.
  friend MultiDiagonal copy(const auto& md)
  {
    MultiDiagonal res;
    for (const auto& [index,offsets,diag] : multidiagonal::range(md)) res.diagonals[index].emplace(offsets,copy(diag));
    return res;
  }

  Diagonals diagonals;

  /// @brief Applies the operator as \f$H/i\f$: accumulates `dpsidt(idxLHS) += diag(idxDiag) * psi(idxRHS)`
  /// for each diagonal.
  void operator () (double, MultiArrayConstView<dcomp,RANK> psi, MultiArrayView<dcomp,RANK> dpsidt) const
  {
    if (diagonals.empty()) return;

#ifndef NDEBUG
    if (psi.extents != calculateAndCheckDimensions(*this)) throw std::runtime_error("Mismatch between StateVector and MultiDiagonal dimensions");
#endif // NDEBUG

    for (const auto& [index,offsets,diag] : multidiagonal::range(this)) {

      Dimensions ubound{psi.extents}, lboundLHS, lboundRHS;
      for (auto&& [i,o,u,lL,lR] : std::views::zip(std::views::iota(0uz,RANK),offsets,ubound,lboundLHS,lboundRHS)) {
        u -= index[i] ? o : 0; lL = index[i] ? 0 : o; lR = index[i] ? o : 0;
      }

      const auto increment=[&] (auto n, const auto& inc, Dimensions& idxLHS, Dimensions& idxRHS, Dimensions& idxDiag)
      {
        using namespace hana::literals;
        if constexpr (n!=0_c) {
          if (idxLHS[n]==ubound[n]-1) {idxLHS[n]=lboundLHS[n]; idxRHS[n]=lboundRHS[n]; idxDiag[n]=0; inc(n-1_c,inc,idxLHS,idxRHS,idxDiag);}
          else {idxLHS[n]++; idxRHS[n]++; idxDiag[n]++;}
        }
        else {idxLHS[0]++; idxRHS[0]++; idxDiag[0]++;}
      };

      for ( Dimensions idxLHS{lboundLHS}, idxRHS{lboundRHS}, idxDiag{}; idxLHS[0]!=ubound[0] ; increment( hana::llong_c<RANK-1>, increment, idxLHS, idxRHS, idxDiag ) )
        dpsidt(idxLHS) += diag(idxDiag)*psi(idxRHS);
    }
  }


  /// @brief Tensor (direct) product: produces a `MultiDiagonal<RANK+RANK2>` with concatenated indices and offsets.
  template <size_t RANK2>
  friend auto operator* (const MultiDiagonal<RANK>& md1, const MultiDiagonal<RANK2>& md2)
  {
    using ResultType = MultiDiagonal<RANK+RANK2>;
    ResultType res;
    for (const auto& [index1,diagToIndex1] : md1.diagonals) for (const auto& [index2,diagToIndex2] : md2.diagonals) {
      typename ResultType::Index resIndex{index2.to_string()+index1.to_string()};
      for (const auto& [offsets1,diag1] : diagToIndex1) for (const auto& [offsets2,diag2] : diagToIndex2)
        res.diagonals[resIndex].emplace(concatenate(offsets1,offsets2),directProduct(diag1,diag2));
    }
    return res;
  }


  /// @name Hermitian conjugation
  //@{

  /// @brief In-place Hermitian conjugation.
  /// Flips `Index` bits for non-zero offsets (transpose) and conjugates all diagonal elements.
  MultiDiagonal& hermitianConjugate()
  {
    Diagonals newDiagonals;
    for (auto&& [index,offsets,diag] : multidiagonal::range(this)) {
      Index newIndex{index}; for (size_t i=0; i<RANK; ++i) if (offsets[i]) newIndex.flip(i);
      conj(diag);
      newDiagonals[newIndex].emplace(offsets,std::move(diag));
    }
    diagonals.swap(newDiagonals);
    return *this;
  }

  MultiDiagonal& dagger() {return hermitianConjugate();}

  friend MultiDiagonal hermitianConjugateOf(const MultiDiagonal& md) {MultiDiagonal res(copy(md)); res.hermitianConjugate(); return res;}

  friend MultiDiagonal twoTimesRealPartOf(const MultiDiagonal& md) {return md+hermitianConjugateOf(md);}
  friend MultiDiagonal twoTimesImagPartOf(const MultiDiagonal& md) {return md-hermitianConjugateOf(md);}
  //@}


  /// @name Algebra
  /// @note Boost.Operator cannot be used since `MultiDiagonal` is not copyable.
  //@{
  MultiDiagonal operator-() const
  {
    MultiDiagonal res{copy(this)};
    for (auto&& [index,offsets,diag] : multidiagonal::range(res)) for (dcomp& v : diag.dataStorage()) v*=-1;
    return res;
  }

  MultiDiagonal operator+() const {return copy(this);}

  /// @brief In-place addition. New diagonals are inserted; existing ones accumulate element-wise.
  MultiDiagonal& operator+=(const MultiDiagonal& md)
  {
#ifndef NDEBUG
    if (auto dim{calculateAndCheckDimensions(*this)}, mdDim{calculateAndCheckDimensions(md)};
      dim!=Dimensions{} && mdDim!=Dimensions{} && dim!=mdDim)
        throw std::runtime_error("Dimension mismatch in addition of MultiDiagonals");
#endif // NDEBUG
    for (const auto& [index,offsets,diag] : multidiagonal::range(md))
      if (auto insertResult=diagonals[index].emplace(offsets,copy(diag)); !insertResult.second)
        for (auto&& [to,from] : std::views::zip(insertResult.first->second.dataStorage(),diag.dataStorage())) to+=from;
    return *this;
  }

  MultiDiagonal& operator-=(const MultiDiagonal& md) {return operator+=(-md);}

  MultiDiagonal& operator*=(scalar auto d)
  {
    for (auto&& [index,offsets,diag] : multidiagonal::range(this)) for (dcomp& v : diag.dataStorage()) v*=d;
    return *this;
  }

  MultiDiagonal& operator/=(scalar auto d) {(*this)*=1./d; return *this;}

  friend auto operator+(const MultiDiagonal& md1, const MultiDiagonal& md2) {MultiDiagonal res{copy(md1)}; res+=md2; return res;}
  friend auto operator-(const MultiDiagonal& md1, const MultiDiagonal& md2) {MultiDiagonal res{copy(md1)}; res-=md2; return res;}

  friend auto operator*(scalar auto v, const MultiDiagonal& md) {MultiDiagonal res{copy(md)}; res*=v; return res;}
  friend auto operator/(const MultiDiagonal& md, scalar auto v) {MultiDiagonal res{copy(md)}; res/=v; return res;}
  //@}


  /// @brief Derives and validates the Hilbert space dimensions from diagonal extents and offsets.
  /// Returns `Dimensions{}` for an empty operator. Throws if any diagonal is inconsistent.
  friend Dimensions calculateAndCheckDimensions(const MultiDiagonal& md)
  {
    if (md.diagonals.empty()) return {};

    auto transformExtentsOfDiagonals = [] (const auto& diagonal) -> Dimensions {
      auto res{diagonal.second.extents};
      for (auto&& [v,o] : std::views::zip(res,diagonal.first)) v+=o;
      return res;
    };

    auto res{transformExtentsOfDiagonals(*md.diagonals.cbegin()->second.cbegin())};

    if (!std::ranges::fold_left_first( std::views::join(md.diagonals | std::views::values) | std::views::transform([=] (const auto& diagonal) {
      return std::ranges::equal(res,transformExtentsOfDiagonals(diagonal));
    }), std::logical_and{} ).value_or(true)) throw std::runtime_error("Dimensions mismatch in MultiDiagonal");

    return res;
  }

  /// @brief Boost.JSON serialization: emits a JSON object keyed by direction-pattern bitstring.
  friend void tag_invoke( const json::value_from_tag&, json::value& jv, const MultiDiagonal& md )
  {
    for (const auto& d : md.diagonals) jv.as_object().emplace( d.first.to_string(), json::value_from(d.second) );
  }

};


/// @brief Composition of two rank-1 `MultiDiagonal` operators: \f$(A|B)_{nm} = \sum_k A_{nk} B_{km}\f$.
///
/// Result diagonals appear at all pairwise offset sums; see `MultiDiagonal.cc` for the four upper/lower cases.
///
/// @note Only rank-1 composition is currently implemented.
/// @note `operator*` has higher precedence than `operator|` in C++; use parentheses in mixed expressions.
MultiDiagonal<1> operator|(const MultiDiagonal<1>&, const MultiDiagonal<1>&);


namespace multidiagonal {

/// @brief Constructs the rank-1 identity operator of dimension @p dim.
MultiDiagonal<1> identity(size_t dim);

} // multidiagonal


} // quantumoperator