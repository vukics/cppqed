// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_QUANTUMOPERATOR_TRIDIAGONAL_H_INCLUDED
#define CPPQEDCORE_QUANTUMOPERATOR_TRIDIAGONAL_H_INCLUDED


#ifndef   QUANTUMOPERATOR_TRIDIAGONAL_MAX_RANK
#define   QUANTUMOPERATOR_TRIDIAGONAL_MAX_RANK BLITZ_ARRAY_LARGEST_RANK
#endif // QUANTUMOPERATOR_TRIDIAGONAL_MAX_RANK

#include "DimensionsBookkeeper.h"
#include "LazyDensityOperator.h"
#include "Types.h"

#include "BlitzTinyOfArrays.tcc"
#include "Operators.h"



/// Comprises modules representing operators of special structure (tridiagonal, sparse) over Hilbert spaces of arbitrary arity
namespace quantumoperator {



//////////////
//
// Tridiagonal
//
//////////////


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
template<int RANK> 
class Tridiagonal 
  : public DimensionsBookkeeper<RANK>, 
    private linalg::VectorSpace<Tridiagonal<RANK> >
{
public:
  static const int LENGTH=tmptools::Power_v<3,RANK>; ///< The number of Diagonals the class has to store
  static const int N_RANK=RANK;                      ///< Reports the class’s rank for example towards DirectProduct.

  typedef blitzplusplus::TinyOfArrays<dcomp,RANK,LENGTH> Diagonals; ///< The class is implemented in terms of a blitzplusplus::TinyOfArrays, this is the class used to store the Diagonals

  typedef typename Diagonals::T_numtype Diagonal; ///< A unary complex blitz array

private:
  typedef DimensionsBookkeeper<RANK> Base;
  
public:
  typedef typename Base::Dimensions Dimensions; ///< Inherited from DimensionsBookkeeper

  typedef quantumdata::StateVectorLow<RANK> StateVectorLow;
  
private:
  typedef blitz::TinyVector<blitz::TinyVector<blitz::Range,3>,RANK> Ranges;
  
  /// \name Conveniencies for member-function signature definitions
  //@{
  static const Diagonal empty;//=Diagonal    ();
  //@}
  
public:
  /// \name Constructors, assignment
  //@{
    /// Constructor for unary Tridiagonal
    /**
     * This is the principal way to create an object of this class, which can be used in the unary case only, as ensured by the trailing dummy argument.
     * This creates an object corresponding to the elementary operator
     * \f[H^\text{elem}/i=\sum_{n=0}^{N-1}\alpha^{(0)}_n\ket{n}\bra{n}+\sum_{n=0}^{N-1-K}\lp\alpha^{(-)}_n\ket{n+K}\bra{n}+\alpha^{(+)}_n\ket{n}\bra{n+K}\rp\f]
     * when `toFreq=false` and
     * \f[H\Int^\text{elem}(t=0)/i=\sum_{n=0}^{N-1-K}\lp e^{\delta_{n}t}\alpha^{(-)}_n\ket{n+K}\bra{n}+e^{-\delta_{n}t}\alpha^{(+)}_n\ket{n}\bra{n+K}\rp\f]
     * when `toFreq=true`. In this case, the \f$\delta\f$ frequencies are calculated out of \f$\alpha^{(0)}\f$ as \f$\delta_{n}=\alpha^{(0)}_{n}-\alpha^{(0)}_{n+K}\f$.
     * 
     * Either of the three initializing arrays might be zero-size, which signifies that the corresponding diagonal is zero, however, if they are nonzero-size,
     * then their sizes must be compatible with each other. If both `minus` and `plus` are zero-size (purely diagonal matrix), then `k` might be zero as well.
     * Violation is detected at runtime, and an exception of type TridiagonalConsistencyErrorException is thrown.
     * 
     * \note It is a dilemma whether the parameter `k` should be considered a compile-time or a runtime parameter.
     * In the majority of cases it is known already an compile time (e.g. ladder operators, angular momentum operators, etc.).
     * The reason why it is treated as a runtime parameter is spatial degrees of freedom. There, operators like \f$sin(kx)\f$, \f$cos(kx)\f$, etc.,
     * are also tridiagonal in momentum space, and we wanted to have the possibility of specifying \f$k\f$ at runtime.
     *
     */
  template<int R=RANK>
  explicit Tridiagonal(std::enable_if_t<R==1,const Diagonal&> zero=empty, ///< corresponds to \f$\alpha^{(0)}\f$
                       size_t k=0,                  ///< corresponds to \f$K\f$ above
                       const Diagonal& minus=empty, ///< corresponds to \f$\alpha^{(-)}\f$
                       const Diagonal& plus=empty,  ///< corresponds to \f$\alpha^{(+)}\f$
                       bool toFreqs=false);         ///< governs whether the main diagonal (\f$\alpha^{(0)}\f$) is converted to frequencies \f$\delta_{n}\f$

    /// Copy constructor with deep-copy semantics
  Tridiagonal(const Tridiagonal& tridiag) : Tridiagonal(tridiag,tridiag.diagonals_,tridiag.differences_,tridiag.tCurrent_,tridiag.freqs_) {}

    /// Constructs the object as the direct product of two Tridiagonal matrices.
    /**
     * This is rather non-trivial, and the calculation of the resulting diagonals’ position is at the moment calculated at runtime,
     * though it could partly be done at compile time. The eventual frequencies are also composed, direct product translates to a “direct sum” in this case.
     * This really makes sense only if the time instants of the two tridiagonals are the same. Violation is detected at runtime,
     * and an exception of type TridiagonalTimeMismatchException is thrown.
     * 
     * All tridiagonals of `RANK>1` in the framework originate from direct products of unary tridiagonals.
     * 
     * \tparam RANK2 the arity of one of the arguments (the other being `RANK-RANK2`)
     */
  template<int RANK2>
  Tridiagonal(const Tridiagonal<RANK2>&, const Tridiagonal<RANK-RANK2>&);
  //@}

  /// \name Class-specific functionality
  //@{
    /// Furnishes a unary tridiagonal with frequencies calculated from `mainDiagonal`.
    /**
     * Note that the matrix may have its own main diagonal in this case, which remains untouched.
     * (This is actually exploited if we want to transform out only part of a main diagonal with interaction picture.)
     * 
     * \return `*this`
     * 
     * \note Works in the unary case only
     * 
     * \todo Extend to arbitrary arity
     */
  template <int R=RANK>
  std::enable_if_t<R==1,Tridiagonal&> furnishWithFreqs(const Diagonal& mainDiagonal);
  
    /// “Applies” the tridiagonal matrix on the state vector `psiprime`, in the vein of structure::Hamiltonian::addContribution(), that is \f$\ket{\Psi'}+=T\ket\Psi\f$.
    /**
     * Finding out which of the 3 to the power of arity diagonals corresponds to which state-vector slice when “applying”, is a formidable task for higher arity,
     * and a significant portion of this calculation is done at compile time. The structure of this problem naturally maps to a recursion. There are 2 possibilities:
     * -# Explicit specializations of the function apply(). In this case, explicitly specialized code must be generated with the Python script 
     * `quantumoperator/applyImpl/generate.py` for each `RANK` (the Python script takes the rank as its first parameter in the command line).
     * The user is encouraged to try it out to see the structure of the problem. Here, the recursion inherent in the problem is shifted to this Python script.
     * This solution sidesteps template metaprogramming, however, it has the drawback that the amount of code such produced grows exponentially with the arity.
     * This possibility is used if the #DO_CONSIDER_EXPLICITLY_SPECIALIZED_TRIDIAGONAL_APPLIES macro **is defined** @ compilation.
     * -# The function is implemented using the recursive doApply, whose second version is present to break the recursivity. In this case, only this second version of the doApply
     * function must be explicitly specialized, which results in much less code.
     * This possibility is used if the #DO_CONSIDER_EXPLICITLY_SPECIALIZED_TRIDIAGONAL_APPLIES macro **is not defined** @ compilation.
     * 
     * \note (Perhaps somewhat surprisingly,) compile-time resource requirement is larger in the 1st case.
     */
  void apply(const StateVectorLow& psi, StateVectorLow& dpsidt) const;

    /// Updates the elements of the matrix to time instant `t` using the stored frequencies.
    /**
     * \return `*this`
     */
  Tridiagonal& propagate(double t);
  
    /// Calculates the quantum average of a Tridiagonal in a quantum state described by a (unary) quantumdata::LazyDensityOperator
    /**
     * \note Not implemented.
     * 
     * \todo Implement.
     */
  template<int R=RANK>
  std::enable_if_t<R==1,dcomp> average(const quantumdata::LazyDensityOperator<1>&, const typename quantumdata::LazyDensityOperator<1>::Idx&);
  //@}

  /// \name Getters
  //@{
  const Diagonals & get           () const {return   diagonals_;} ///< returns the diagonals
  const Dimensions& getDifferences() const {return differences_;}

  double            getTime       () const {return    tCurrent_;}

  const Diagonals & getFreqs      () const {return       freqs_;}
  //@}


  /// \name Algebra
  //@{
  Tridiagonal& operator*=(dcomp dc); ///< Naively implemented, could be templated if need arises – Frequencies untouched.
  Tridiagonal& operator/=(dcomp dc) {(*this)*=1./dc; return *this;} ///< ”
  Tridiagonal& operator*=(double d) {(*this)*=dcomp(d,0); return *this;} ///< ”
  Tridiagonal& operator/=(double d) {(*this)*=1./dcomp(d,0); return *this;} ///< ”

  // void           hermitianConjugateSelf();
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
  const Tridiagonal hermitianConjugate    () const;
  const Tridiagonal dagger                () const {return hermitianConjugate();} ///< Same as hermitianConjugate

    /// Naive addition.
    /**
    * The structures of the two objects must match, and there is a rather complicated set of rules as to what is meant by “matching”:
    * - Tridiagonal freely mixes with a purely diagonal matrix throught addition.
    * - If both have off-diagonals as well, the `k`s must match.
    * - If both store frequencies, they must be equal.
    * 
    * Any violation of these rules is detected at runtime, and an exception of type TridiagonalStructureMismatchException is thrown.
    */
  Tridiagonal& operator+=(const Tridiagonal&);

  const Tridiagonal operator-() const {return Tridiagonal(*this,blitzplusplus::negate(diagonals_),differences_,tCurrent_,freqs_);} ///< Returns a (deep) copy with negated diagonals. Frequencies remain untouched.

  const Tridiagonal operator+() const {return *this;} ///< Returns a (deep) copy.
  
  Tridiagonal& operator-=(const Tridiagonal& tridiag) {(*this)+=-tridiag; return *this;} ///< Implemented in terms of operator+=
  //@}
  
private:
  Tridiagonal(const Base& base, const Diagonals& diagonals, const Dimensions& differences, double tCurrent, const Diagonals& freqs)
    : Base(base), diagonals_(blitzplusplus::DeepCopy(),diagonals), differences_(differences), tCurrent_(tCurrent), freqs_(blitzplusplus::DeepCopy(),freqs) {}

  /// \name Apply implementation
  //@{
  template<int START, typename V_DPSIDT, typename V_A, typename V_PSI, int REMAINING>
  void doApply(tmptools::integral_c<REMAINING>,const Ranges&, const StateVectorLow&, StateVectorLow&) const; ///< Recursive implementation.

  template<int START, typename V_DPSIDT, typename V_A, typename V_PSI>
  void doApply(tmptools::integral_c<0>,const Ranges&, const StateVectorLow&, StateVectorLow&) const; ///< “Specialization” for the `REMAINING=0` case to break the recursion: This should be specialized for all `RANK`s in Tridiagonal.tcc. We rely on the Boost.Preprocessor library, in particular the \refBoostConstruct{BOOST_PP_ITERATE,preprocessor/doc/ref/iterate.html} macro.

  struct FillRangesHelper
  {
    typedef const typename StateVectorLow::T_index Bound;
    
    FillRangesHelper(Ranges& ranges, const Bound& ubound, const Dimensions& k) : ranges_(ranges), ubound_(ubound), k_(k) {}

    template<typename ICW> void operator()(ICW);
  
  private:
    Ranges& ranges_;
    const Bound& ubound_;
    const Dimensions& k_;

  };

  const Ranges fillRanges(const typename StateVectorLow::T_index&) const; ///< Helper to apply
  //@}

  /// \name Data
  //@{
  Diagonals diagonals_;

  Dimensions differences_;

  double tCurrent_;

  Diagonals freqs_;
  //@}

};


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


} // quantumoperator


namespace structure { namespace freesystem { using Tridiagonal=quantumoperator::Tridiagonal<1>; } }


#endif // CPPQEDCORE_QUANTUMOPERATOR_TRIDIAGONAL_H_INCLUDED
