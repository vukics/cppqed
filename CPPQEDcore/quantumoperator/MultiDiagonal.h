// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "StateVector.h"

#include <bitset>


/// Comprises modules representing operators of special structure (multidiagonal, sparse) over Hilbert spaces of arbitrary arity
namespace quantumoperator {

using namespace ::quantumdata;

namespace multidiagonal {


auto range(auto&& md)
{
  if constexpr (std::is_pointer_v<std::decay_t<decltype(md)>>) return range(*md);
  else return md.diagonals | std::views::transform([&](auto&& outer_entry) {
    return outer_entry.second | std::views::transform([&](auto&& inner_entry) {
      return std::tie(outer_entry.first, inner_entry.first, inner_entry.second);
    });
  }) | std::views::join;
}

template <typename MD>
void for_each(MD&& md, auto&& func)
{
  for (auto&& [index,offsets,diag] : range(std::forward<MD>(md)))
    func(std::forward<decltype(index)>(index),
         std::forward<decltype(offsets)>(offsets),
         std::forward<decltype(diag)>(diag));
}


} // multidiagonal


/// A sparse matrix storing a set of diagonals
/**
 * non-trivial functionality
 * - application as Hamiltonian (after propagation) on StateVectorView
 * - addition (new diagonals may appear!)
 * - composition
 * - direct product
 * - (composition with Sigma)
 *
 * TODO: MultiDiagonal should be applicable also as a jump operator, but there are serious difficulties in indexing.
 *
 * Diagonal class :
 * - as multiarray:
 *   * factors
 *   * optionally frequencies
 * - optionally as a single lambda: an envelope function
 *
 * + a global envelope can also be useful in certain cases
 *
 * The envelopes can be evaluated when applying the object as a Hamiltonian
 *
 */
template <size_t RANK>
struct MultiDiagonal
{
  using Dimensions = Extents<RANK>;
  using Index = std::bitset<RANK>; // true means upper (or the main) diagonal, false lower
  using Offsets = Dimensions;
  using Diagonal = MultiArray<dcomp,RANK>;
  using DiagToIdx = std::map<Offsets,Diagonal>;

  /// \note `std::bitset` has hash but no comparison, whereas `std::array` has lexicographic comparison by default
  using Diagonals = std::unordered_map<Index,DiagToIdx>;

  MultiDiagonal(const MultiDiagonal&) = delete; MultiDiagonal& operator=(const MultiDiagonal&) = delete;
  MultiDiagonal(MultiDiagonal&&) = default; MultiDiagonal& operator=(MultiDiagonal&&) = default;

  explicit MultiDiagonal(auto&&... args) : diagonals{std::forward<decltype(args)>...} {}

  friend MultiDiagonal copy(auto&& md)
  {
    MultiDiagonal res;
    for (const auto& [index,offsets,diag] : multidiagonal::range(std::forward<decltype(md)>(md))) res.diagonals[index].emplace(offsets,copy(diag));
    return res;
  }

  Diagonals diagonals;
  // double tCurrent=0;

  /// Applying as a Hamiltonian
  void operator () (double t, StateVectorConstView<RANK> psi, StateVectorView<RANK> dpsidt) const
  {
    if (diagonals.empty()) return;

#ifndef NDEBUG
    if (psi.extents != calculateAndCheckDimensions(*this)) throw std::runtime_error("Mismatch between StateVector and MultiDiagonal dimensions");
#endif // NDEBUG

    for (const auto& [index,offsets,diag] : multidiagonal::range(this))  {

      Dimensions ubound{psi.extents}, lboundLHS, lboundRHS;
      for (auto&& [i,o,u,lL,lR] : std::views::zip(std::views::iota(0uz,RANK),offsets,ubound,lboundLHS,lboundRHS)) {
        u -= index[i] ? o : 0; lL = index[i] ? 0 : o; lR = index[i] ? o : 0;
      }

      // Here, a generalized version of `cppqedutils::incrementMultiIndex` could be used, but let’s not touch that at the moment
      const auto increment=[&] (auto n, const auto& inc, Dimensions& idxLHS, Dimensions& idxRHS, Dimensions& idxDiag)
      {
        using namespace hana::literals;
        if constexpr (n!=0_c) {
          if (idxLHS[n]==ubound[n]-1) {idxLHS[n]=lboundLHS[n]; idxRHS[n]=lboundRHS[n]; idxDiag[n]=0; inc(n-1_c,inc,idxLHS,idxRHS,idxDiag);}
          else {idxLHS[n]++; idxRHS[n]++; idxDiag[n]++;}
        }
        else {idxLHS[0]++; idxRHS[0]++; idxDiag[0]++;} // This will eventually put the index into an illegal state, but this is how all (unchecked) iterators work.
      };

      for ( Dimensions idxLHS{lboundLHS}, idxRHS{lboundRHS}, idxDiag{}; idxLHS[0]!=ubound[0] ; increment( hana::llong_c<RANK-1>, increment, idxLHS, idxRHS, idxDiag ) )
        dpsidt(idxLHS) += diag(idxDiag)*psi(idxRHS) ;
    }
  }


  /// Direct product
  template <size_t RANK2>
  friend auto operator* (const MultiDiagonal<RANK>& md1, const MultiDiagonal<RANK2>& md2)
  {
    using ResultType = MultiDiagonal<RANK+RANK2>;
    ResultType res;
    for (const auto& [index1,diagToIndex1] : md1.diagonals) for (const auto& [index2,diagToIndex2] : md2.diagonals) {
      typename ResultType::Index resIndex{index1.to_string()+index2.to_string()};
      for (const auto& [offsets1,diag1] : diagToIndex1) for (const auto& [offsets2,diag2] : diagToIndex2)
        res.diagonals[resIndex].emplace(concatenate(offsets1,offsets2),directProduct(diag1,diag2));
    }
    return res;
  }


  /// \name Hermitian conjugation
  //@{
  MultiDiagonal& hermitianConjugate()
  {
    Diagonals newDiagonals;
    for (auto&& [index,offsets,diag] : multidiagonal::range(this)) {
        /// in the original index all those elements must be negated for which offset is nonzero
        Index newIndex{index}; for (size_t i=0; i<RANK; ++i) if (offsets[i]) newIndex.flip(i);
        conj(diag); // conjugate the diagonal
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

  /// calculates a propagator from `mainDiagonal` and stores them in `frequencies`
  // MultiDiagonal& furnishWithFreqs(const Diagonal& mainDiagonal) requires (RANK==1);
  
  /// propagates the elements of the matrix to time instant `t` using `frequencies`
  // MultiDiagonal& propagate(double t);
  
  /// Calculates the quantum average of a MultiDiagonal in a quantum state described by a (unary) quantumdata::LazyDensityOperator
  // template<int R=RANK> std::enable_if_t<R==1,dcomp> average(const quantumdata::LazyDensityOperator<1>&, const typename quantumdata::LazyDensityOperator<1>::Idx&);


  /// \name Algebra
  /** We cannot use old techniques like Boost.Operator, since MultiDiagonal is not copyable */
  //@{
  MultiDiagonal operator-() const
  {
    MultiDiagonal res{copy(this)};
    for (auto&& [index,offsets,diag] : multidiagonal::range(res)) for (dcomp& v : diag.dataStorage()) v*=-1;
    return res;
  }

  MultiDiagonal operator+() const {return copy(this);}

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

  friend void to_json( json& jv, const MultiDiagonal& md )
  {
    for (const auto& d : md.diagonals) jv.emplace(d.first.to_string(),d.second);
  }

};



/// Composition
/**
 * It’s unfortunate that operator* has precedence over operator| according to c++ precedence rules,
 * so parentheses must be used in expressions containing both of them, but we can live with this.
 */
MultiDiagonal<1> operator|(const MultiDiagonal<1>&, const MultiDiagonal<1>&);


namespace multidiagonal {

MultiDiagonal<1> identity(size_t dim);

} // multidiagonal

/*
/// A free-standing version of Tridiagonal::apply \related Tridiagonal
template<int RANK>
inline void
apply(const quantumdata::StateVectorLow<RANK>& psi, quantumdata::StateVectorLow<RANK>& dpsidt, const Tridiagonal<RANK>& tridiag)
{
  tridiag.apply(psi,dpsidt);
}


/// Same as Tridiagonal::furnishWithFreqs, but returns a copy of its first argument furnished with frequencies \related Tridiagonal
template<int RANK>
Tridiagonal<RANK>
furnishWithFreqs(const Tridiagonal<RANK>& tridiag,                        ///< Tridiagonal whose copy is to be furnished with frequencies
                 const typename Tridiagonal<RANK>::Diagonal& mainDiagonal ///< main diagonal determining the frequencies to be used for furnishing
                );


/// Unary zero operator as a Tridiagonal \related Tridiagonal
const Tridiagonal<1> zero    (size_t);
/// Unary identity operator as a Tridiagonal \related Tridiagonal
const Tridiagonal<1> identity(size_t);



template<int RANK>
const typename Tridiagonal<RANK>::Diagonal Tridiagonal<RANK>::empty;



/// Direct product \related Tridiagonal
template<int RANK1, int RANK2>
inline
Tridiagonal<RANK1+RANK2>
operator*(const Tridiagonal<RANK1>& t1, const Tridiagonal<RANK2>& t2)
{
  return Tridiagonal<RANK1+RANK2>(t1,t2);
}


/// Returns the anti-Hermitian operator \f$T-T^\dagger\f$ \related Tridiagonal
template<int RANK>
inline 
Tridiagonal<RANK>
tridiagMinusHC     (const Tridiagonal<RANK>& tridiag) {return tridiag-tridiag.hermitianConjugate();}


/// Returns the Hermitian operator \f$T+T^\dagger\f$ \related Tridiagonal
template<int RANK>
inline 
Tridiagonal<RANK>
tridiagPlusHC      (const Tridiagonal<RANK>& tridiag) {return tridiag+tridiag.hermitianConjugate();}


/// Returns the anti-Hermitian operator \f$(T+T^\dagger)/i\f$ \related Tridiagonal
template<int RANK>
inline 
Tridiagonal<RANK>
tridiagPlusHC_overI(const Tridiagonal<RANK>& tridiag) {return tridiagPlusHC(tridiag)/1i;}


/// \related Tridiagonal
template<int RANK>
std::ostream& operator<<(std::ostream&, const Tridiagonal<RANK>&);
*/

} // quantumoperator


/// Class representing a (unary) tridiagonal matrix or direct product of such matrices
/**
 * Let us consider the following situation, which displays the full potential of the class: We have a system consisting of \f$M\f$ subsystems,
 * each featuring a free (in general, non-Hermitian) Hamiltonian, which is diagonal in the given system’s Hilbert space.
 * In addition, we have one more term of a special form in the Hamiltonian coupling all the subsystems (we are considering \f$H/i\f$,
 * because this is what appears in the framework as noted @ structure::Hamiltonian::addContribution):
 * \f[H/i=\lp H^\text{free}+H^\text{interaction}\rp/i=\sum_{m=0}^{M-1}\bigotimes_{k=0}^{m-1}\mathbf{1}_k\otimes\lp\sum_{n_m=0}^{N_m-1}\omega_{m,n}\ket{n_m}\bra{n_m}\rp\otimes\bigotimes_{k=m+1}^{M-1}\mathbf{1}_k+\bigotimes_{m=0}^{M-1}\lp\sum_{n_m=0}^{N_m-1}\alpha^{(0)}_{m,n}\ket{n_m}\bra{n_m}\right.+\left.\sum_{n_m=0}^{N_m-1-K_m}\lp\alpha^{(+)}_{m,n}\ket{n_m}\bra{n_m+K_m}+\alpha^{(-)}_{m,n}\ket{n_m+K_m}\bra{n_m}\rp\rp.\f]
 * Here, the coefficients \f$\omega\f$ and \f$\alpha\f$ are in general complex with the dimension of frequency (\f$\hbar=1\f$).
 * The \f$N_m\f$s are the dimensions of the subsystems’ Hilbert spaces in which the vectors \f$\ket{n_m}\f$ form an orthonormal basis.
 *
 * The Hamiltonian is indeed in a special form because the interaction term is a direct product of such operators acting on the Hilbert spaces of the individual subsystems,
 * whose matrix contains only three diagonals of nonzero elements. Hence the name of the class Tridiagonal, although this usually refers to the case when \f$K=1\f$.
 *
 * Now let us transform to the interaction picture defined by \f$H^\text{free}\f$. The Hamiltonian in interaction picture reads
 * \f[H\Int(t)/i=\bigotimes_{m=0}^{M-1}\lp\sum_{n_m=0}^{N_m-1}\alpha^{(0)}_{m,n}\ket{n_m}\bra{n_m}+\sum_{n_m=0}^{N_m-1-K_m}\lp e^{\delta_{m,n}t}\alpha^{(-)}_{m,n}\ket{n_m+K_m}\bra{n_m}+e^{-\delta_{m,n}t}\alpha^{(+)}_{m,n}\ket{n_m}\bra{n_m+K_m}\rp\rp,\f]
 * where \f$\delta_{m,n}=\omega_{m,n}-\omega_{m,n+K_m}\f$.
 *
 * Quite generally, the class is designed to store and manipulate Hamiltonians of the form either \f$H\f$ or \f$H\Int(t)\f$ above for an arbitrary number of subsystems.
 * In the latter case, it also stores and manipulates the \f$\delta_{m,n}\f$ frequencies. In particular, the Hamiltonian can be evaluated at a given time \f$t\f$,
 * applied on a state vector and combined with other Tridiagonal instances using algebraic operations.
 *
 * Tridiagonal internally bookkeeps to which time instant its state corresponds.
 *
 * Policy of storing diagonals: only the non-zeros. Policy of storing frequencies: either none or all of them. Because of these policies,
 * Tridiagonal can represent a diagonal matrix without significant overhead.
 *
 * \note The frequencies \f$\omega_{m,n}\f$ are not necessarily real here, Tridiagonal works also for non-Hermitian operators.
 * On non-unitary transformations in quantum mechanics cf. [these notes](http://optics.szfki.kfki.hu/~vukics/Pictures.pdf).
 *
 * \tparamRANK corresponding to \f$M\f$ above
 *
 * \note A serious limitation of Tridiagonal is that the composition of two such operators does not in general yield an operator of the same form.
 * This is one of the reasons why we are planning to deprecate this class in favour of the much more general form
 * \f[H\Int(t)=\bigotimes_{m=0}^{M-1}\sum_{i_m\in\mathbb{K}_m}\sum_{n_m=0}^{N_m-1-i_m}e^{\delta^{i_m}_{m,n}t}\alpha^{i_m}_{m,n}\ket{n_m+i_m}\bra{n_m},\f]
 * with \f$\mathbb{K}_m=\left\{K_m^{(0)},K_m^{(1)},…\right\}\f$ an arbitrary set, and \f$\delta_{m,n}^{(i_m)}=i\lp\omega_{m,n+i_m}-\omega_{m,n}\rp\f$.
 *
 */
