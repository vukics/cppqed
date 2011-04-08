// -*- C++ -*-
#ifndef   LAZY_DENSITY_OPERATOR_SMART_ITERATOR_IMPL_INCLUDED
#define   LAZY_DENSITY_OPERATOR_SMART_ITERATOR_IMPL_INCLUDED

#include "DensityOperator.h"
#include "StateVector.h"

#include "Algorithm.h"
#include "BlitzArraySliceIterator.h"
#include "Functional.h"


namespace quantumdata {


namespace ldo {


namespace details {


template<int RANK, typename V> struct ExtendV : mpl::fold<V,
							  V,
							  mpl::push_back<mpl::_1,
									 mpl::plus<mpl::_2,
										   mpl::int_<RANK> > > >
{};


class OutOfRange           : public cpputils::Exception {};
class NoSuchImplementation : public cpputils::Exception {};


template<typename V>
class DI_ImplSpecial
{
public:
  static const int RANK=mpl::size<V>::value;

  typedef LazyDensityOperator<RANK>    LDO         ;
  typedef boost::shared_ptr<const LDO> LDO_SmartPtr;


  DI_ImplSpecial(const LDO&    , mpl:: true_) : ldoSmartPtr_()             , isEnd_( true) {}
  // In this case, the SmartPtr is never actually touched

  DI_ImplSpecial(const LDO& ldo, mpl::false_) : ldoSmartPtr_(dispatch(ldo)), isEnd_(false) {}


  const LDO_SmartPtr dispatch(const LDO& ldo)
  {
    using blitzplusplus::basi::Transposer;

    typedef StateVector    <RANK> SV;
    typedef DensityOperator<RANK> DO;

    if      (const SV*const sV=dynamic_cast<const SV*>(&ldo)) {
      typename SV::StateVectorLow temp(sV->operator()()); 
      // We do not want to transpose that StateVectorLow which is the storage of sV.
      return LDO_SmartPtr(new SV(Transposer<RANK,V>::transpose(temp),byReference));
    }
    else if (const DO*const dO=dynamic_cast<const DO*>(&ldo)) {
      typename DO::DensityOperatorLow temp(dO->operator()());
      return LDO_SmartPtr(new DO(Transposer<2*RANK,typename ExtendV<RANK,V>::type>::transpose(temp),byReference));
    }
    else throw NoSuchImplementation();
  }

  bool isEqualTo(const DI_ImplSpecial& other) const {return isEnd_==other.isEnd_;}

  void doIncrement() {if (!isEnd_) isEnd_=true; else throw OutOfRange();}

  const LDO& dereference() const {if (isEnd_) throw OutOfRange(); return *ldoSmartPtr_;}

  virtual ~DI_ImplSpecial() {}

private:
  mutable LDO_SmartPtr ldoSmartPtr_;

  bool isEnd_;

};



template<int RANK, typename V>
class DI_Impl
{
public:
  typedef blitzplusplus::basi::Iterator<RANK,V,true> BASI;
  typedef typename BASI::Impl MII;

  typedef TTD_LAZY_DENSITY_OPERATOR_RES LazyDensityOperatorRes;

  bool isEqualTo(const DI_Impl& other) const {return getMII()==other.getMII();}

  virtual void doIncrement() = 0;
  virtual const LazyDensityOperatorRes& dereference() const = 0;

  virtual ~DI_Impl() {}

private:
  virtual const MII& getMII() const = 0;

};



////////////////
//
// DI_SV_Impl
//
////////////////



template<int RANK, typename V>
class DI_SV_Impl : public DI_Impl<RANK,V>
{
public:
  typedef DI_Impl<RANK,V> Base;

  typedef typename Base::BASI Impl;
  typedef typename Base::MII  MII ;

  typedef StateVector<mpl::size<V>::value>        StateVectorRes        ;
  typedef boost::shared_ptr<const StateVectorRes> StateVectorResSmartPtr;

  typedef typename Base::LazyDensityOperatorRes LazyDensityOperatorRes;

  template<bool IS_END>
  DI_SV_Impl(const StateVector<RANK>& psi, mpl::bool_<IS_END> tag) : impl_(psi(),tag), stateVectorResSmartPtr_() {}

private:
  void doIncrement() {++impl_;}

  const MII& getMII() const {return impl_();}

  const LazyDensityOperatorRes& dereference() const {stateVectorResSmartPtr_.reset(new StateVectorRes(*impl_,byReference)); return *stateVectorResSmartPtr_;}

  Impl impl_;

  mutable StateVectorResSmartPtr stateVectorResSmartPtr_;

};


////////////////
//
// DI_DO_Impl
//
////////////////



template<int RANK, typename V>
class DI_DO_Impl : public DI_Impl<RANK,V>
{
public:
  typedef DI_Impl<RANK,V> Base;

  typedef typename Base::MII MII;

  static const int RANKRES=mpl::size<V>::value;

  typedef DensityOperator<RANKRES>                     DensityOperatorRes        ;
  typedef boost::shared_ptr<const DensityOperatorRes > DensityOperatorResSmartPtr;


  typedef typename Base::LazyDensityOperatorRes LazyDensityOperatorRes;

  typedef typename ExtendV<RANK,V>::type ExtendedV;
  
  typedef blitzplusplus::basi::Iterator<2*RANK,ExtendedV,true> BASI;

  typedef typename BASI::CA    DensityOperatorLow   ;
  typedef typename BASI::CARes DensityOperatorLowRes;

  typedef typename MII::IdxTiny IdxTinyHalf;
  typedef typename TTD_IDXTINY(2*IdxTinyHalf::numElements) IdxTiny;

  DI_DO_Impl(const DensityOperator<RANK>& rho, mpl::false_);
  DI_DO_Impl(const DensityOperator<RANK>& rho, mpl:: true_);

private:
  template<bool TAG>
  static MII ctorHelper(const DensityOperator<RANK>&);

  void doIncrement() {++mii_;}

  const MII& getMII() const {return mii_;}

  const LazyDensityOperatorRes& dereference() const;

  MII mii_;

  mutable DensityOperatorLow    densityOperatorLow_   ; 
  mutable DensityOperatorLowRes densityOperatorLowRes_;

  mutable DensityOperatorResSmartPtr densityOperatorResSmartPtr_; 

};


template<int RANK, typename V>
template<bool TAG>
typename DI_DO_Impl<RANK,V>::MII
DI_DO_Impl<RANK,V>::ctorHelper(const DensityOperator<RANK>& rho)
{
  using blitzplusplus::basi::details::FilterOut;
  using blitzplusplus::halfCutTiny;
  return MII(FilterOut<RANK,V>(halfCutTiny(rho().lbound())),
	     FilterOut<RANK,V>(halfCutTiny(rho().ubound())),
	     mpl::bool_<TAG>());
}


template<int RANK, typename V>
DI_DO_Impl<RANK,V>::DI_DO_Impl(const DensityOperator<RANK>& rho, mpl::false_)
  : mii_(ctorHelper<false>(rho)),
    densityOperatorLow_(), densityOperatorLowRes_(), densityOperatorResSmartPtr_()
{
  densityOperatorLow_.reference(rho());
  BASI::transpose(densityOperatorLow_);
}


template<int RANK, typename V>
DI_DO_Impl<RANK,V>::DI_DO_Impl(const DensityOperator<RANK>& rho, mpl:: true_)
  : mii_(ctorHelper< true>(rho)),
    densityOperatorLow_(), densityOperatorLowRes_(), densityOperatorResSmartPtr_()    
{
}

template<int RANK, typename V>
const typename DI_DO_Impl<RANK,V>::LazyDensityOperatorRes&
DI_DO_Impl<RANK,V>::dereference() const
{
  using namespace blitzplusplus;
  IdxTinyHalf idxHalf(*mii_);
  IdxTiny idx(concatenateTinies(idxHalf,idxHalf));

  densityOperatorResSmartPtr_.reset(new DensityOperatorRes(
							   BASI::index(
								       densityOperatorLow_,
								       densityOperatorLowRes_,
								       idx
								       ),
							   byReference
							   )
				    );
  return *densityOperatorResSmartPtr_;
}


//////////////////
//
// Implementations
//
//////////////////


template<int RANK, typename V, bool IS_END>
typename DiagonalIterator<RANK,V>::Impl
makeDI_Impl(const LazyDensityOperator<RANK>& ldo, mpl::false_)
// The last argument governing whether the special implementation is needed.
{
  static const mpl::bool_<IS_END> isEnd=mpl::bool_<IS_END>();
  typedef typename DiagonalIterator<RANK,V>::Impl DI_ImplSmartPtr;

  if      (const StateVector    <RANK>*const stateVector    =dynamic_cast<const StateVector    <RANK>*>(&ldo))
    return DI_ImplSmartPtr(new DI_SV_Impl<RANK,V>(*stateVector    ,isEnd));
  else if (const DensityOperator<RANK>*const densityOperator=dynamic_cast<const DensityOperator<RANK>*>(&ldo))
    return DI_ImplSmartPtr(new DI_DO_Impl<RANK,V>(*densityOperator,isEnd));
  else
    throw NoSuchImplementation();
  
}


template<int RANK, typename V, bool IS_END>
typename DiagonalIterator<RANK,V>::Impl
makeDI_Impl(const LazyDensityOperator<RANK>& ldo, mpl::true_)
// The last argument governing whether the special implementation is needed.
{
  static const mpl::bool_<IS_END> isEnd=mpl::bool_<IS_END>();
  typedef typename DiagonalIterator<RANK,V>::Impl DI_ImplSmartPtr;

  return DI_ImplSmartPtr(new DI_ImplSpecial<V>(ldo,isEnd));

}


} // details


template<int RANK, typename V>
template<bool IS_END>
DiagonalIterator<RANK,V>::DiagonalIterator(const LazyDensityOperator<RANK>& ldo, mpl::bool_<IS_END>)
  : impl_(details::makeDI_Impl<RANK,V,IS_END>(ldo,mpl::bool_<IS_SPECIAL>()))
{}



template<int RANK, typename V>
DiagonalIterator<RANK,V>&
DiagonalIterator<RANK,V>::operator++()
{
  impl_->doIncrement();
  return *this;
}

template<int RANK, typename V>
const typename DiagonalIterator<RANK,V>::LazyDensityOperatorRes&
DiagonalIterator<RANK,V>::operator*() const
{
  return impl_->dereference();
}


template<int RANK, typename V>
bool
DiagonalIterator<RANK,V>::isEqualTo(const DiagonalIterator& other) const
{
  return impl_->isEqualTo(*other.impl_);
}


} // ldo


template<int RANK> template<typename V>
const ldo::DiagonalIterator<RANK,V>
LazyDensityOperator<RANK>::begin(V) const
{
  return ldo::DiagonalIterator<RANK,V>(*this,mpl::false_());
}

template<int RANK> template<typename V>
const ldo::DiagonalIterator<RANK,V>
LazyDensityOperator<RANK>::end  (V) const
{
  return ldo::DiagonalIterator<RANK,V>(*this,mpl:: true_());
}


template<int RANK, typename F, typename V, typename T>
const T
partialTrace(const LazyDensityOperator<RANK>& matrix, F function, V v, T)
{
  ldo::DiagonalIterator<RANK,V> begin(matrix.begin(v));
  T init(function(*begin));

  using namespace cpputils;
  return T(
	   accumulate(++begin,
		      matrix.end(v),
		      init,
		      function,
		      cpputils::plus<T>()
		      )
	   );
}


} // quantumdata



#endif // LAZY_DENSITY_OPERATOR_SMART_ITERATOR_IMPL_INCLUDED
