// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_QUANTUMDATA_LAZYDENSITYOPERATOR_H_INCLUDED
#define CPPQEDCORE_QUANTUMDATA_LAZYDENSITYOPERATOR_H_INCLUDED

#include "ArrayBase.h" // for ByReference
#include "QuantumDataFwd.h"

#include "DimensionsBookkeeper.h"

#include "BlitzArray.h"
#include "BlitzTinyExtensions.h"
#include "ComplexExtensions.h"
#include "MathExtensions.h"
#include "SliceIterator.tcc"
#include "MultiIndexIterator.h"

#include <boost/operators.hpp>

#include <boost/mpl/size.hpp>
#include <boost/mpl/if.hpp>

#include <boost/range/numeric.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <memory>

namespace mpl=boost::mpl;


/// Comprises classes representing the state of composite quantum systems and providing various interfaces to manipulate this data
/** Some of its most important classes fit into a single class rooted in the virtual interface LazyDensityOperator. */
namespace quantumdata {


/// Common interface for calculating quantum averages
/**
 * In a quantum-simulation framework, users should be able to write code for calculating quantum expectation values from quantumdata, 
 * independently of whether this data is represented by state vectors or density operators, in orthogonal or non-orthogonal bases. 
 * One obvious solution is relying on the formula
 * \f[\avr{A}=\Tr{A\rho}\f]
 * (\f$A\f$ being an observable and \f$\rho\f$ the density operator of the system), to write code only for the density-operator case,
 * and fall back to this in the state-vector case as well, by calculating a dyad from the state vector. This is, however, extremely wasteful, 
 * since usually not all the matrix elements of \f$\rho\f$ are needed for calculating the average, furthermore, for large dimensionality,
 * this solution may become outright unaffordable in terms of memory: for large systems, we may afford to store \f$\ket\Psi\f$, but not \f$\ket\Psi\bra\Psi\f$.
 * 
 * The solution adopted for this problem in the framework is represented by LazyDensityOperator, which provides a common interface 
 * for all the four cases StateVector, DensityOperator, and their non-orthogonal counterparts, to calculate quantum averages from their data.
 * The “laziness” amounts to that in the case of state vectors, only those elements of the density operator are calculated that are actually asked for.
 * 
 * This class is totally immutable, all its member functions being constant.
 * 
 * \tparamRANK
 */
template<int RANK> 
class LazyDensityOperator 
  : public DimensionsBookkeeper<RANK,true>,
    std::enable_shared_from_this<const LazyDensityOperator<RANK>>
{
public:
  typedef std::shared_ptr<const LazyDensityOperator> Ptr; ///< Many class templates in the framework define shared pointers to their own types, in a template-metafunction like manner
  
  typedef DimensionsBookkeeper<RANK,true> Base;

  typedef typename Base::Dimensions Dimensions; ///< Inherited from DimensionsBookkeeper

  typedef IdxTiny<RANK> Idx; ///< The type used for indexing the “rows” and the “columns”: a tiny vector of integers (multi-index)

  virtual ~LazyDensityOperator() {}

private:
  class IndexerProxy
  {
  public:
    IndexerProxy(const LazyDensityOperator* ldo, const Idx& firstIndex) : ldo_(ldo), firstIndex_(firstIndex) {}

    const dcomp operator()(const Idx& secondIndex) const {return ldo_->index(firstIndex_,secondIndex);}

    template<typename... SubscriptPack>
    const dcomp operator()(int s0, SubscriptPack... subscriptPack) const
    {
      static_assert( mpl::size<mpl::vector<SubscriptPack...> >::value==RANK-1 , "Incorrect number of subscripts for LazyDensityOperator::IndexerProxy." );
      return operator()(Idx(s0,subscriptPack...));
    }

    operator double() const {return real(ldo_->index(firstIndex_,firstIndex_));}

  private:
    const LazyDensityOperator*const ldo_;
    const Idx firstIndex_;

  };

public:
  /// Multi-matrix style indexing via Idx type
  const IndexerProxy operator()(const Idx& firstIndex) const {return IndexerProxy(this,firstIndex);}

  /// Multi-matrix style indexing via packs of integers
  /**
   * This allows for the very convenient indexing syntax (e.g. a ternary LazyDensityOperator `matrix` indexed by multi-indeces `i` and `j`):
   *
   *     matrix(i0,i1,i2)(j0,j1,j2)
   *
   * while
   *
   *     matrix(i0,i1,i2)
   *
   * returns a proxy implicitly convertible to a `double`, giving the diagonal element corresponding to the multi-index `i`
   *
   * \note The number of indeces in the multi-index is checked @ compile time.
   *
   */
  template<typename... SubscriptPack>
  const auto operator()(int s0, SubscriptPack... subscriptPack) const
  {
    static_assert( sizeof...(SubscriptPack)==RANK-1 , "Incorrect number of subscripts for LazyDensityOperator." );
    return operator()(Idx(s0,subscriptPack...));
  }

  double trace() const {return trace_v();} ///< Returns the trace (redirected to a pure virtual)
  
  /// \name Slicing-related functionality
  //@{
    /// Return the ldo::DiagonalIterator corresponding to the beginning of the sequence of slices defined by `V`
    /** Cf. \ref slicinganldo "rationale"
     * \tparamV
     */
  template<typename V>
  inline auto begin() const
  {
    return DiagonalIterator<V>(*this,mpl::false_());
  }

  template<typename V>
  inline auto end  () const ///< ” for the end
  {
    return DiagonalIterator<V>(*this,mpl:: true_());
  }
  //@}
  
protected:
  LazyDensityOperator(const Dimensions& dims) : Base(dims) {}

private:
  virtual const dcomp index(const Idx& i, const Idx& j) const = 0;
  
  virtual double trace_v() const = 0;
  
  /// Iterator for slices of a LazyDensityOperator that are diagonal in the dummy indices
  /**
  * Cf. \ref slicinganldo "rationale"
  * 
  * \tparamV
  * 
  * It's inherently a const iterator since LazyDensityOperator is immutable.
  * Models an [InputIterator](http://www.cplusplus.com/reference/std/iterator/InputIterator/), implemented with the help of
  * \refBoostConstruct{input_iterator_helper,utility/operators.htm#iterator} from Boost.Operator.
  * 
  */
  template<typename V>
  class DiagonalIterator;
 

};


/// The primary tool for performing \ref slicinganldo "slice iterations"
/**
 * On higher levels of the framework (cf. eg. BinarySystem, Composite), this function is used exclusively for performing LazyDensityOperator slice iteration.
 * 
 * \tparamV
 * \tparam T an arithmetic type which must be default-constructible
 * \tparam F a callable type with signature `const T(const typename ldo::DiagonalIterator<RANK,V>::value_type&)`. This signature is equivalent to `const T(const LazyDensityOperator<mpl::size<V>::value>&)`.
 * 
 */
template<typename V, typename T, int RANK, typename F>
const T
partialTrace(const LazyDensityOperator<RANK>& matrix, F function)
{
  auto begin(matrix.template begin<V>());
  T init(function(*begin)); // we take a way around default constructing T, so that the implicit interface is less stringently defined

  using namespace boost;
  return accumulate(make_iterator_range(++begin,matrix.template end<V>()) | adaptors::transformed(function) , init );
}


/// Turns the data of a LazyDensityOperator into a real 1D array
/**
 * The result is an array starting with the diagonals (real) of the matrix, optionally followed by the off-diagonals (real & imaginary) of the upper triangle.
 * For this to make sense, the matrix is assumed to be Hermitian.
 * 
 * \param offDiagonals governs whether the upper triangle of the (multi-)matrix is included in the result
 */
template<int RANK>
const DArray<1> deflate(const LazyDensityOperator<RANK>& matrix, bool offDiagonals)
{
  using mathutils::sqr;
  
  const size_t dim=matrix.getTotalDimension();
  
  DArray<1> res(offDiagonals ? sqr(dim) : dim);
  
  typedef cpputils::MultiIndexIterator<RANK> Iterator;
  const Iterator etalon(matrix.getDimensions()-1,cpputils::mii::begin);
  
  size_t idx=0;

  for (Iterator i(etalon); idx<dim; ++i)
    res(idx++)=matrix(*i);
  
  if (offDiagonals)
    for (Iterator i(etalon); idx<sqr(dim); ++i)
      for (Iterator j=++Iterator(i); j!=etalon.getEnd(); ++j) {
        dcomp matrixElement(matrix(*i)(*j));
        res(idx++)=real(matrixElement);
        res(idx++)=imag(matrixElement);
      }

  return res;

}


/** \class LazyDensityOperator
 * 
 * \Semantics
 * 
 * - *Unary system*: Assume a mode represented in Fock basis with ladder-operator \f$a\f$. To calculate the quantum expectation value
 * \f[\avr{a^2}=\Tr{a^2\rho}=\sum_i\sqrt{i(i-1)}\,\rho_{i;i-2},\f] one can write the following function:
 *   ~~~
 *   include "LazyDensityOperator.h"
 *   
 *   const dcomp calculateASqr(const LazyDensityOperator<1>& matrix)
 *   {
 *     dcomp res;
 *     for (int i=2; i<matrix.getTotalDimension(); ++i) res+=sqrt(i*(i-1))*matrix(i,i-2);
 *     return res;
 *   }
 *   ~~~
 * 
 * - *Binary system*: Assume two modes represented in Fock bases with ladder-operators \f$a\f$ and \f$b\f$, respectively.
 *   To calculate the quantum expectation value\f[\avr{a^\dag b}=\Tr{a^\dag b\rho}=\sum_{i,j}\sqrt{(i+1)j}\,\rho_{(i,j);(i+1,j-1)},\f]
 *   one can write the following function:
 *   ~~~ 
 *   include "LazyDensityOperator.h"
 * 
 *   const dcomp calculateADaggerB(const LazyDensityOperator<2>& matrix)
 *   {
 *     const LazyDensityOperator<2>::Dimensions dim(matrix.getDimensions());
 * 
 *     dcomp res;
 *     for (int i=0; i<dim[0]-1; ++i) for (int j=1; j<dim[1]; ++j) res+=sqrt((i+1)*j)*matrix(i,j)(i+1,j-1);
 *     return res;
 *   }
 *   ~~~
 * 
 */

/** \fn template<typename V, typename T, int RANK, typename F> const T partialTrace(const LazyDensityOperator<RANK>&, F function)
 * 
 * \Semantics
 * 
 * The function iterates through all the combinations of the dummy indices, for each slice it takes the value returned by the functor `function`, and accumulates these values.
 * In the following we give some instructive examples of usage.
 * 
 * - *Calculating the full partial density operator of a unary subsystem*: 
 *   The partial density operator of any unary subsystem of a system of any rank may be calculated as follows (e.g. from a StateVector):
 *   ~~~
 *   template<int SUBSYSTEM, int RANK> // for the subsystem indexed by the index SUBSYSTEM
 *   const DensityOperator<1>
 *   partialTraceOfUnarySubsystem(const StateVector<RANK>& psi)
 *   {
 *     return partialTrace<tmptools::Vector<SUBSYSTEM>,DensityOperator<1> >(psi,densityOperatorize<1>);
 *   }
 *   ~~~
 * 
 * - *Calculating correlations for two harmonic-oscillator modes*: Assume having the function `calculateADaggerB` as defined @ LazyDensityOperator.
 *   Then, if these modes are embedded at index positions 3 and 1 (note that the order might be important) in a larger system of arity larger than 3,
 *   we can write the following function:
 *   ~~~
 *   template<int RANK> // RANK>3
 *   const dcomp
 *   calculateADaggerB_atPositions3and1(const LazyDensityOperator<RANK>& matrix)
 *   {
 *     return partialTrace<tmptools::Vector<3,1>,dcomp>(matrix,calculateADaggerB);
 *   }
 *   ~~~
 *   This will calculate \f$\avr{a^\dag b}\f$ (now \f$b\f$ and \f$a\f$ being the ladder operators of the modes at position 1 and 3, respectively)
 *   for the partial density operator of the embedded system (but without explicitly calculating this partial density operator).
 *
 * - *Accumulating an arbitrary set of quantum averages*: An arbitrary set of real or complex numbers can be stored in a DArray<1> or CArray<1>
 *   (or even `std::valarray`s), as these types fulfill all the requirements on type `T`. The former was used to define the abstract interface structure::Averaged
 *   one of its implementations being :class:`Mode`. E.g. a function for a harmonic-oscillator mode calculating \f$\avr{a^\dag a}\f$, \f$\avr{\lp a^\dag a\rp^2}\f$,
 *   and the real and imaginary parts of \f$\avr{a}\f$, may be defined as
 *   ~~~
 *   typedef DArray<1> Averages;
 *
 *   const Averages
 *   calculateModeAverages(const LazyDensityOperator<1>& matrix)
 *   {
 *     Averages averages(4); averages=0;
 *
 *     for (int n=1; n<int(matrix.getDimension()); n++) {
 *       double diag=matrix(n);
 *       averages(0)+=  n*diag;
 *       averages(1)+=n*n*diag;
 *
 *       double sqrtn=sqrt(n);
 *       dcomp offdiag(matrix(n,n-1));
 *       averages(2)+=sqrtn*real(offdiag);
 *       averages(3)+=sqrtn*imag(offdiag);
 *     }
 *
 *     return averages;
 *
 *   }
 *   ~~~
 *   Then the following function will calculate these averages for a mode embedded in a larger system at index position `MODE_POSITION`:
 *   ~~~
 *   template<int RANK, int MODE_POSITION> // for the subsystem indexed by the index MODE_POSITION
 *   const Averages
 *   calculateEmbeddedModeAverages(const StateVector<RANK>& psi)
 *   {
 *     return partialTrace(psi,calculateModeAverages,tmptools::Vector<MODE_POSITION>(),Averages()); 
 *   }
 *   ~~~
 */


} // quantumdata

/** \page slicinganldo Slicing a LazyDensityOperator
 * 
 * \tableofcontents
 * 
 * Analogously to \ref multiarrayconcept "slicing state vectors", it is also necessary to slice LazyDensityOperator objects
 * because for calculating quantum expectation values of subsystem-observables (e.g. in Composite objects),
 * the partial-trace density operator is needed.
 * 
 * For the partial trace, however, only such elements of the full density operator are needed as are diagonal in the indices *not* belonging
 * to the given subsystem (dummy indices). This is the reason why the tool performing the iteration is called DiagonalIterator.
 *
 * Slicing is fully recursive, a sliced LazyDensityOperator (usually obtained by dereferencing a DiagonalIterator) can be further sliced.
 * 
 * Notes on implementation {#slicinganldoimplementation}
 * =======================
 * 
 * Difficulty: LazyDensityOperator is an abstract interface, and the organization of its data (and hence the actual procedure of slicing)
 * varies along implementations.
 * 
 * \image html ldoDiagonalIterator.png
 * 
 * Implemented using a classical inheritance-based strategy idiom, together with both compile-time and run-time implementation selection
 * (similarly to blitzplusplus::basi::Iterator, a special implementation (DiagonalIterator::DI_ImplSpecial) is needed when the size
 * of the compile-time vector `V` equals `RANK`).
 * 
 * The run-time polymorphy of LazyDensityOperator necessitates run-time implementation selection.
 * This happens when the actual implementation (`DI_SV_Impl` or `DI_DO_Impl`, or their non-orthogonal counterparts) is chosen 
 * depending on whether the given LazyDensityOperator is a StateVector or a DensityOperator (or their non-orthogonal counterparts).
 * 
 */


#define INPUT_IteratorHelper boost::input_iterator_helper<DiagonalIterator<V>,const LazyDensityOperator<mpl::size<V>::value> >

template<int RANK> template<typename V>
class quantumdata::LazyDensityOperator<RANK>::DiagonalIterator
  : public INPUT_IteratorHelper 
{
private:
  typedef INPUT_IteratorHelper Base;

#undef INPUT_IteratorHelper

  typedef typename Base::value_type LazyDensityOperatorRes; ///< Equivalent to `const LazyDensityOperator<mpl::size<V>::value>`

public:
  static const bool IS_SPECIAL=(RANK==mpl::size<V>::value); ///< Signifies whether the special implementation is needed

  /// Constructor
  /**
  * Similarly to blitzplusplus::basi::Iterator, it can be initialised either to the beginning or to the end of the sequence.
  * \tparam IS_END governs the end-ness
  */
  template<bool IS_END>
  DiagonalIterator(const LazyDensityOperator& ldo, mpl::bool_<IS_END>) : impl_(makeImpl<IS_END>(ldo,mpl::bool_<IS_SPECIAL>())) {}
  
  /// \name Necessary members of an input iterator
  //@{
  DiagonalIterator& operator++() {impl_->doIncrement(); return *this;} ///< Immediately delegated to the implementation

  const LazyDensityOperatorRes& operator*() const {return impl_->dereference();} ///< ”
  
  bool operator==(const DiagonalIterator& other) const {return impl_->isEqualTo(*other.impl_);} ///< ”
  //@}
  
  class NoSuchImplementation : public cpputils::Exception {};
  
  class DI_ImplSpecial ///< Special implementation in the case `RANK=mpl::size<V>::value`
  {
  public:
    typedef LazyDensityOperator LDO;
    typedef std::shared_ptr<const LDO> LDO_Ptr;

    class OutOfRange : public cpputils::Exception {};

    DI_ImplSpecial(const LDO&    , cpputils::mii::End  ) : ldoPtr_()             , isEnd_( true) {}
    // In this case, the Ptr is never actually touched

    DI_ImplSpecial(const LDO& ldo, cpputils::mii::Begin) : ldoPtr_(dispatch(ldo)), isEnd_(false) {}


    const LDO_Ptr dispatch(const LDO& ldo)
    {
      using blitzplusplus::basi::Transposer;

      typedef StateVector    <RANK> SV;
      typedef DensityOperator<RANK> DO;

      if      (const auto sV=dynamic_cast<const SV*>(&ldo)) {
        typename SV::StateVectorLow temp(sV->getArray());
        // We do not want to transpose that StateVectorLow which is the storage of sV.
        return std::make_shared<SV>(Transposer<RANK,V>::transpose(temp),byReference);
      }
      else if (const auto dO=dynamic_cast<const DO*>(&ldo)) {
        typename DO::DensityOperatorLow temp(dO->getArray());
        return std::make_shared<DO>(Transposer<2*RANK,typename tmptools::ExtendVector<RANK,V>::type>::transpose(temp),byReference);
      }
      else throw NoSuchImplementation();
    }

    bool isEqualTo(const DI_ImplSpecial& other) const {return isEnd_==other.isEnd_;}

    void doIncrement() {if (!isEnd_) isEnd_=true; else throw OutOfRange();}

    const LDO& dereference() const {if (isEnd_) throw OutOfRange(); return *ldoPtr_;}

  private:
    mutable LDO_Ptr ldoPtr_;

    bool isEnd_;

  };

  class DI_Impl ///< Base class for non-special implementations
  {
  public:
    typedef blitzplusplus::basi::Iterator<RANK,V,true> BASI;
    typedef typename BASI::Impl MII;

    typedef typename DiagonalIterator::LazyDensityOperatorRes LazyDensityOperatorRes;

    bool isEqualTo(const DI_Impl& other) const {return getMII()==other.getMII();}

    virtual void doIncrement() = 0;
    virtual const LazyDensityOperatorRes& dereference() const = 0;

    virtual ~DI_Impl() {}

  private:
    virtual const MII& getMII() const = 0;

  };

  /// Pointer to implementation
  /**
  * We are using the classical inheritance-based pointer-to-implementation technique, together with some compile-time dispatching.
  * \refBoostConstruct{eval_if,mpl/doc/refmanual/eval-if.html} from Boost.MPL here guarantees that `DI_ImplSpecial` gets instantiated
  * only in the special case
  */
  typedef std::shared_ptr<typename mpl::eval_if_c<IS_SPECIAL,
                                                  mpl::identity<DI_ImplSpecial>,
                                                  mpl::identity<DI_Impl       > 
                                                  >::type
                          > Impl;

private:
  Impl impl_;

  class DI_SV_Impl : public DI_Impl
  {
  public:
    typedef typename DI_Impl::BASI BASI;
    typedef typename DI_Impl::MII MII;

    typedef StateVector<mpl::size<V>::value> StateVectorRes;

    template<bool IS_END>
    DI_SV_Impl(const StateVector<RANK>& psi, mpl::bool_<IS_END> tag) : basi_(psi.getArray(),tag), stateVectorResPtr_() {}

  private:
    void doIncrement() {++basi_;}

    const MII& getMII() const {return basi_();}

    const LazyDensityOperatorRes& dereference() const {stateVectorResPtr_.reset(new StateVectorRes(*basi_,byReference)); return *stateVectorResPtr_;}

    BASI basi_;

    mutable std::shared_ptr<const StateVectorRes> stateVectorResPtr_;

  };

  
  class DI_DO_Impl : public DI_Impl
  {
  public:
    typedef typename DI_Impl::MII MII;
    
    static const int RANKRES=mpl::size<V>::value;

    typedef DensityOperator<RANKRES> DensityOperatorRes;

    typedef typename tmptools::ExtendVector<RANK,V>::type ExtendedV;
    
    typedef blitzplusplus::basi::Iterator<2*RANK,ExtendedV,true> BASI;

    typedef typename BASI::CA    DensityOperatorLow   ;
    typedef typename BASI::CARes DensityOperatorLowRes;
 
    DI_DO_Impl(const DensityOperator<RANK>& rho, cpputils::mii::Begin)
      : mii_(ctorHelper<false>(rho)), densityOperatorLow_(), densityOperatorLowRes_(), densityOperatorResPtr_()
    {
      densityOperatorLow_.reference(rho.getArray());
      BASI::transpose(densityOperatorLow_);
    }

    DI_DO_Impl(const DensityOperator<RANK>& rho, cpputils::mii::End  )
      : mii_(ctorHelper< true>(rho)), densityOperatorLow_(), densityOperatorLowRes_(), densityOperatorResPtr_()    
    {}
    
  private:
    template<bool IS_END>
    static MII ctorHelper(const DensityOperator<RANK>& rho)
    {
      using blitzplusplus::basi::filterOut;
      using blitzplusplus::halfCutTiny;
      
      return MII(filterOut<RANK,V>(halfCutTiny(rho.getArray().lbound())),
                 filterOut<RANK,V>(halfCutTiny(rho.getArray().ubound())),
                 mpl::bool_<IS_END>());
    }

    void doIncrement() {++mii_;}

    const MII& getMII() const {return mii_;}

    const LazyDensityOperatorRes& dereference() const
    {
      return *(densityOperatorResPtr_=std::make_shared<DensityOperatorRes>(BASI::index(densityOperatorLow_,
                                                                                       densityOperatorLowRes_,
                                                                                       blitzplusplus::concatenateTinies(*mii_,*mii_)),
                                                                           byReference));
    }

    MII mii_;

    mutable DensityOperatorLow    densityOperatorLow_   ; 
    mutable DensityOperatorLowRes densityOperatorLowRes_;

    mutable std::shared_ptr<const DensityOperatorRes> densityOperatorResPtr_; 

  };

  
  template <bool IS_END>
  static auto makeImpl(const LazyDensityOperator& ldo, mpl::true_)
  {
    return std::make_shared<DI_ImplSpecial>(ldo,mpl::bool_<IS_END>());
  }
  
  template <bool IS_END>
  static std::shared_ptr<DI_Impl> makeImpl(const LazyDensityOperator& ldo, mpl::false_)
  {
    static const auto isEnd{mpl::bool_<IS_END>()};

    if      (const auto stateVector    =dynamic_cast<const StateVector    <RANK>*>(&ldo))
      return std::make_shared<DI_SV_Impl>(*stateVector    ,isEnd);
    else if (const auto densityOperator=dynamic_cast<const DensityOperator<RANK>*>(&ldo))
      return std::make_shared<DI_DO_Impl>(*densityOperator,isEnd);
    else throw NoSuchImplementation();
  }
  
};




#endif // CPPQEDCORE_QUANTUMDATA_LAZYDENSITYOPERATOR_H_INCLUDED
