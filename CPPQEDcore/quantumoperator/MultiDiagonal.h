// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "StateVector.h"



/// Comprises modules representing operators of special structure (multidiagonal, sparse) over Hilbert spaces of arbitrary arity
namespace quantumoperator {


/// A sparse matrix storing a set of diagonals
/**
 * non-trivial functionality
 * - application as Hamiltonian (after propagation) on StateVectorView
 * - addition (new diagonals may appear!)
 * - composition
 * - direct product
 * - (composition with Sigma)
 *
 * \todo
 */
template <size_t RANK>
struct MultiDiagonal
{
  using Dimensions=::cppqedutils::Extents<RANK>;

  using Offsets = std::array<int,RANK>;

  struct CompareOffsets
  {
    bool operator() (const Offsets& o1, const Offsets o2) const {
      auto recurse = [&] (Offsets::const_iterator i1, Offsets::const_iterator i2, const auto r) {
        if (i1==o1.end()) return false;
        else if (*i1==*i2) return r(++i1,++i2,r);
        else return *i1<*i2;
      };
      return recurse(o1.begin(),o2.begin(),recurse);
    }
  };

  using Diagonal = ::cppqedutils::MultiArray<dcomp,RANK>;

  /// indexing via an indexer ensures that all diagonals are indexed 0 … e-1-o, where e is the extent along the given dimension, and o is the offset
  template <typename D> requires std::same_as<std::decay_t<D>,Diagonal>
  static auto indexer(const Offsets& offsets, D&& diagonal) {
    return [&] (::cppqedutils::Extents<RANK> idx) {
      for (auto&& [o,i] : std::views::zip(offsets,idx)) if (o<0) i+=o;
      return diagonal(idx);
    };
  }

  // need to add comparison for std::arrays
  using Diagonals = std::map<Offsets,Diagonal,CompareOffsets>;

  MultiDiagonal(const MultiDiagonal&) = delete; MultiDiagonal& operator=(const MultiDiagonal&) = delete;
  MultiDiagonal(MultiDiagonal&&) = default; MultiDiagonal& operator=(MultiDiagonal&&) = default;

  explicit MultiDiagonal(Dimensions dimensions, auto&&... args) : dimensions{dimensions}, diagonals{std::forward<decltype(args)>...} {}

  Dimensions dimensions;
  Diagonals diagonals;
  // double tCurrent=0;

  void hermitianConjugateSelf() {
    Diagonals newDiagonals;
    for (auto iter=diagonals.begin(); iter!=diagonals.end(); ++iter) {
      Offsets newOffsets{};
      for (size_t i=0; i<RANK; ++i) newOffsets[i]=-iter->first[i]; // negate the offsets
      conj(iter->second); // conjugate the diagonal
      newDiagonals.insert(std::pair{newOffsets,std::move(iter->second)});
    }
    diagonals.swap(newDiagonals);
  }

  /// calculates a propagator from `mainDiagonal` and stores them in `frequencies`
  // MultiDiagonal& furnishWithFreqs(const Diagonal& mainDiagonal) requires (RANK==1);
  
  /// propagates the elements of the matrix to time instant `t` using `frequencies`
  // MultiDiagonal& propagate(double t);
  
  /// Calculates the quantum average of a MultiDiagonal in a quantum state described by a (unary) quantumdata::LazyDensityOperator
  // template<int R=RANK> std::enable_if_t<R==1,dcomp> average(const quantumdata::LazyDensityOperator<1>&, const typename quantumdata::LazyDensityOperator<1>::Idx&);


  /// \name Algebra
  //@{
  MultiDiagonal& operator*=(dcomp dc); ///< Naively implemented, could be templated if need arises – Frequencies untouched.
  MultiDiagonal& operator/=(dcomp dc) {(*this)*=1./dc; return *this;} ///< ”
  MultiDiagonal& operator*=(double d) {(*this)*=dcomp(d,0); return *this;} ///< ”
  MultiDiagonal& operator/=(double d) {(*this)*=1./dcomp(d,0); return *this;} ///< ”

    /// Returns a newly constructed object, which is the Hermitian conjugate of `this`.
    /**
    * Transposition involves a non-trivial permutation of diagonals, which could be done @ compile-time, but at the moment it’s runtime.
    * 
    * Frequencies need not be transposed, because they are identical to their transpose.
    * 
    * \note Eventually, an in-place version of Hermitian conjugation could be also implemented, if need for this arises.
    * 
    * \todo Implement an in-place version.
    */
/*  const MultiDiagonal hermitianConjugate    () const;
  const MultiDiagonal dagger                () const {return hermitianConjugate();} ///< Same as hermitianConjugate

  MultiDiagonal& operator+=(const MultiDiagonal&);

  const MultiDiagonal operator-() const {return MultiDiagonal(*this,blitzplusplus::negate(diagonals_),differences_,tCurrent_,freqs_);} ///< Returns a (deep) copy with negated diagonals. Frequencies remain untouched.

  const MultiDiagonal operator+() const {return *this;} ///< Returns a (deep) copy.
  
  MultiDiagonal& operator-=(const MultiDiagonal& tridiag) {(*this)+=-tridiag; return *this;} ///< Implemented in terms of operator+=
  //@}
 */
};


MultiDiagonal<1> compose(const MultiDiagonal<1>&, const MultiDiagonal<1>);

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
