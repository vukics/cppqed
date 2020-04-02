// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDCORE_QUANTUMDATA_LAZYDENSITYOPERATORSLICEITERATOR_TCC_INCLUDED
#define   CPPQEDCORE_QUANTUMDATA_LAZYDENSITYOPERATORSLICEITERATOR_TCC_INCLUDED

#include "LazyDensityOperatorSliceIterator.h"

#include "DensityOperator.tcc"
#include "LazyDensityOperator.h"
#include "StateVector.h"

#include "BlitzArraySliceIterator.tcc"
#include "BlitzTinyExtensions.tcc"

#include <boost/make_shared.hpp>


namespace quantumdata {


namespace ldo {

using namespace cpputils::mii;

class NoSuchImplementation : public cpputils::Exception {};

namespace details {


template<int RANK, typename V>
struct ExtendV : mpl::fold<V,
                           V,
                           mpl::push_back<mpl::_1,
                                          mpl::plus<mpl::_2,
                                                    mpl::int_<RANK>
                                                    >
                                          >
                           >
{};


} // details


template<int RANK, typename V>
class DiagonalIterator<RANK,V>::DI_ImplSpecial
{
public:
  typedef LazyDensityOperator<RANK>    LDO    ;
  typedef boost::shared_ptr<const LDO> LDO_Ptr;

  class OutOfRange : public cpputils::Exception {};

  DI_ImplSpecial(const LDO&    , End  ) : ldoPtr_()             , isEnd_( true) {}
  // In this case, the Ptr is never actually touched

  DI_ImplSpecial(const LDO& ldo, Begin) : ldoPtr_(dispatch(ldo)), isEnd_(false) {}


  const LDO_Ptr dispatch(const LDO& ldo)
  {
    using blitzplusplus::basi::Transposer;

    typedef StateVector    <RANK> SV;
    typedef DensityOperator<RANK> DO;

    if      (const auto sV=dynamic_cast<const SV*>(&ldo)) {
      typename SV::StateVectorLow temp(sV->getArray());
      // We do not want to transpose that StateVectorLow which is the storage of sV.
      return boost::make_shared<SV>(Transposer<RANK,V>::transpose(temp),byReference);
    }
    else if (const auto dO=dynamic_cast<const DO*>(&ldo)) {
      typename DO::DensityOperatorLow temp(dO->getArray());
      return boost::make_shared<DO>(Transposer<2*RANK,typename details::ExtendV<RANK,V>::type>::transpose(temp),byReference);
    }
    else throw NoSuchImplementation();
  }

  bool isEqualTo(const DI_ImplSpecial& other) const {return isEnd_==other.isEnd_;}

  void doIncrement() {if (!isEnd_) isEnd_=true; else throw OutOfRange();}

  const LDO& dereference() const {if (isEnd_) throw OutOfRange(); return *ldoPtr_;}

  virtual ~DI_ImplSpecial() {}

private:
  mutable LDO_Ptr ldoPtr_;

  bool isEnd_;

};



template<int RANK, typename V>
class DiagonalIterator<RANK,V>::DI_Impl
{
public:
  typedef blitzplusplus::basi::Iterator<RANK,V,true> BASI;
  typedef typename BASI::Impl MII;

  typedef typename DiagonalIterator<RANK,V>::LazyDensityOperatorRes LazyDensityOperatorRes;

  bool isEqualTo(const DI_Impl& other) const {return getMII()==other.getMII();}

  virtual void doIncrement() = 0;
  virtual const LazyDensityOperatorRes& dereference() const = 0;

  virtual ~DI_Impl() {}

private:
  virtual const MII& getMII() const = 0;

};


//////////////////
//
// Implementations
//
//////////////////


namespace details {
  

template <int RANK, typename V, bool IS_SPECIAL> class MakeImpl;


template <int RANK, typename V> class MakeImpl<RANK,V,false>
{
private:
  class DI_SV_Impl : public DiagonalIterator<RANK,V>::DI_Impl
  {
  public:
    typedef typename DiagonalIterator<RANK,V>::DI_Impl Base;

    typedef typename Base::BASI Impl;
    typedef typename Base::MII  MII ;

    typedef StateVector<mpl::size<V>::value>        StateVectorRes   ;
    typedef boost::shared_ptr<const StateVectorRes> StateVectorResPtr;

    typedef typename Base::LazyDensityOperatorRes LazyDensityOperatorRes;

    template<bool IS_END>
    DI_SV_Impl(const StateVector<RANK>& psi, mpl::bool_<IS_END> tag) : impl_(psi.getArray(),tag), stateVectorResPtr_() {}

  private:
    void doIncrement() {++impl_;}

    const MII& getMII() const {return impl_();}

    const LazyDensityOperatorRes& dereference() const {stateVectorResPtr_.reset(new StateVectorRes(*impl_,byReference)); return *stateVectorResPtr_;}

    Impl impl_;

    mutable StateVectorResPtr stateVectorResPtr_;

  };
  

  class DI_DO_Impl : public DiagonalIterator<RANK,V>::DI_Impl
  {
  public:
    typedef typename DiagonalIterator<RANK,V>::DI_Impl Base;

    typedef typename Base::MII MII;

    static const int RANKRES=mpl::size<V>::value;

    typedef DensityOperator<RANKRES>                     DensityOperatorRes   ;
    typedef boost::shared_ptr<const DensityOperatorRes > DensityOperatorResPtr;


    typedef typename Base::LazyDensityOperatorRes LazyDensityOperatorRes;

    typedef typename details::ExtendV<RANK,V>::type ExtendedV;
    
    typedef blitzplusplus::basi::Iterator<2*RANK,ExtendedV,true> BASI;

    typedef typename BASI::CA    DensityOperatorLow   ;
    typedef typename BASI::CARes DensityOperatorLowRes;

    DI_DO_Impl(const DensityOperator<RANK>& rho, Begin)
      : mii_(ctorHelper<false>(rho)), densityOperatorLow_(), densityOperatorLowRes_(), densityOperatorResPtr_()
    {
      densityOperatorLow_.reference(rho.getArray());
      BASI::transpose(densityOperatorLow_);
    }

    DI_DO_Impl(const DensityOperator<RANK>& rho, End  )
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
      return *(densityOperatorResPtr_=boost::make_shared<DensityOperatorRes>(BASI::index(densityOperatorLow_,
                                                                                         densityOperatorLowRes_,
                                                                                         blitzplusplus::concatenateTinies(*mii_,*mii_)),
                                                                             byReference));
    }

    MII mii_;

    mutable DensityOperatorLow    densityOperatorLow_   ; 
    mutable DensityOperatorLowRes densityOperatorLowRes_;

    mutable DensityOperatorResPtr densityOperatorResPtr_; 

  };

public:  
  template<bool IS_END>
  static typename DiagonalIterator<RANK,V>::Impl doIt(const LazyDensityOperator<RANK>& ldo)
  {
    using boost::make_shared;

    static const mpl::bool_<IS_END> isEnd=mpl::bool_<IS_END>();

    if      (const auto stateVector    =dynamic_cast<const StateVector    <RANK>*>(&ldo))
      return make_shared<DI_SV_Impl>(*stateVector    ,isEnd);
    else if (const auto densityOperator=dynamic_cast<const DensityOperator<RANK>*>(&ldo))
      return make_shared<DI_DO_Impl>(*densityOperator,isEnd);
    else throw NoSuchImplementation();
    
  }
  
};


template <int RANK, typename V> class MakeImpl<RANK,V,true>
{
public:
  template <bool IS_END>
  static typename DiagonalIterator<RANK,V>::Impl doIt(const LazyDensityOperator<RANK>& ldo)
  {
    return boost::make_shared<typename DiagonalIterator<RANK,V>::DI_ImplSpecial>(ldo,mpl::bool_<IS_END>());
  }
  
};  


} // details


template<int RANK, typename V> template<bool IS_END>
DiagonalIterator<RANK,V>::DiagonalIterator(const LazyDensityOperator<RANK>& ldo, mpl::bool_<IS_END>)
  : impl_(details::MakeImpl<RANK,V,IS_SPECIAL>::template doIt<IS_END>(ldo))
{}


} // ldo


} // quantumdata



#endif // CPPQEDCORE_QUANTUMDATA_LAZYDENSITYOPERATORSLICEITERATOR_TCC_INCLUDED
