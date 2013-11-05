// -*- C++ -*-
/// \briefFileDefault
#ifndef QUANTUMDATA_NONORTHOGONALSTATEVECTOR_H_INCLUDED
#define QUANTUMDATA_NONORTHOGONALSTATEVECTOR_H_INCLUDED

#include "NonOrthogonalStateVectorFwd.h"
#include "NonOrthogonalDensityOperator.h"

#include "StateVector.h"
#include "Transformation.h"



namespace quantumdata {

template<int RANK, typename TRAFO>
class NonOrthogonalStateVector:
    public StateVector<RANK>,
    private linalg::VectorSpace<NonOrthogonalStateVector<RANK, TRAFO> >
{
public:
  typedef StateVector<RANK> Base;

  typedef typename Base::Dimensions Dimensions;
  typedef typename Base::StateVectorLow StateVectorLow;
  typedef typename Base::DensityOperatorLow DensityOperatorLow;

  typedef typename Base::Idx Idx;

  typedef TRAFO Trafo;
  typedef transformation::Traits<Trafo> TrafoTraits;

  using Base::operator(); using Base::vectorView; using Base::getTotalDimension;

  explicit NonOrthogonalStateVector(const StateVectorLow& psi, const Trafo& tf, ByReference)
    : Base(psi,byReference), dual_(psi.shape()), transformation_(tf){}

  explicit NonOrthogonalStateVector(const Dimensions& dimensions,
                                    const Trafo& tf, bool init=true)
    : Base(dimensions, init), dual_(dimensions), transformation_(tf){}

  NonOrthogonalStateVector(const NonOrthogonalStateVector& sv)
    : Base(sv), dual_(sv().shape()),
      transformation_(sv.getTrafo()){}

  double norm() const;

  void update() const;

  const StateVectorLow& dual() const {return dual_;};
  const Trafo& getTrafo() const {return transformation_;}

  template<typename OTHER>
  NonOrthogonalStateVector& operator=(const OTHER& other) {operator()()=other; return *this;}

  // naive operations for vector space

  NonOrthogonalStateVector& operator+=(const NonOrthogonalStateVector& psi) {Base::operator+=(psi); return *this;}
  NonOrthogonalStateVector& operator-=(const NonOrthogonalStateVector& psi) {Base::operator-=(psi); return *this;}

  const NonOrthogonalStateVector operator-() const {return NonOrthogonalStateVector(-operator()());}
  const NonOrthogonalStateVector operator+() const {return *this;} // deep copy

  template<typename OTHER>
  NonOrthogonalStateVector& operator*=(const OTHER& dc) {Base::operator*=(dc); return *this;}

  template<typename OTHER>
  NonOrthogonalStateVector& operator/=(const OTHER& dc) {Base::operator/=(dc); return *this;}

  void addTo(NonOrthogonalDensityOperator<RANK>&) const;

private:
  mutable StateVectorLow dual_;

  const Trafo transformation_;

};



template<typename SV1, typename SV2>
struct TensorType : mpl::identity<NonOrthogonalStateVector<SV1::N_RANK+SV2::N_RANK,
                                                           typename transformation::Compose<typename SV1::Trafo,
                                                                                            typename SV2::Trafo>::type> >
{
};


/** \cond */
template<typename SV1, int RANK2>
struct TensorType<SV1,StateVector<RANK2> > 
  : mpl::identity<NonOrthogonalStateVector<SV1::N_RANK+RANK2,
                                           typename transformation::Compose<typename SV1::Trafo,
                                                                            transformation::Identity<RANK2> >::type> >
{
};


template<int RANK1, typename SV2>
struct TensorType<StateVector<RANK1>,SV2> 
  : mpl::identity<NonOrthogonalStateVector<RANK1+SV2::N_RANK,
                                           typename transformation::Compose<transformation::Identity<RANK1>,
                                                                            typename SV2::Trafo>::type> >
{
};
/** \endcond */



#define RETURN_typeKERNEL transformation::Compose<TRAFO1,TRAFO2>
#define RETURN_type NonOrthogonalStateVector<RANK1+RANK2,typename RETURN_typeKERNEL::type>


template<int RANK1, typename TRAFO1, int RANK2, typename TRAFO2>
inline
const RETURN_type
operator*(const NonOrthogonalStateVector<RANK1,TRAFO1>& sv1, const NonOrthogonalStateVector<RANK2,TRAFO2>& sv2)
{
  using namespace blitzplusplus; 
  return RETURN_type(doDirect<dodirect::multiplication>(sv1(),sv2()),RETURN_typeKERNEL::compose(sv1.getTrafo(),sv2.getTrafo()),byReference);
}

#undef RETURN_typeKERNEL

#define RETURN_typeKERNEL transformation::Compose<TRAFO1,transformation::Identity<RANK2> >

template<int RANK1, typename TRAFO1, int RANK2>
inline
const RETURN_type
operator*(const NonOrthogonalStateVector<RANK1,TRAFO1>& sv1, const StateVector<RANK2>& sv2)
{
  using namespace blitzplusplus; 
  return RETURN_type(doDirect<dodirect::multiplication>(sv1(),sv2()),RETURN_typeKERNEL::compose(sv1.getTrafo(),transformation::Identity<RANK2>()),byReference);
}

#undef RETURN_typeKERNEL

#define RETURN_typeKERNEL transformation::Compose<transformation::Identity<RANK1>,TRAFO2>

template<int RANK1, int RANK2, typename TRAFO2>
inline
const RETURN_type
operator*(const StateVector<RANK1>& sv1, const NonOrthogonalStateVector<RANK2,TRAFO2>& sv2)
{
  using namespace blitzplusplus; 
  return RETURN_type(doDirect<dodirect::multiplication>(sv1(),sv2()),RETURN_typeKERNEL::compose(transformation::Identity<RANK1>(),sv2.getTrafo()),byReference);
}

#undef RETURN_type
#undef RETURN_typeKERNEL


} // quantumdata

#endif // QUANTUMDATA_NONORTHOGONALSTATEVECTOR_H_INCLUDED
