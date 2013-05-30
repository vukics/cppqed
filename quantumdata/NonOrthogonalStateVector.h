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


/** \page nonorthogonalformalism Nonorthogonal basis sets

 * ### Formalism
 * 
 * We adopt the so-called covariant-contravariant formalism (see also \cite{artacho91}), which in physics is primarily known from the theory of relativity.
 * Assume we have a basis \f$\lbr\ket{i}\rbr_{i\in\mathbb{N}}\f$, where the basis vectors are nonorthogonal, so that the metric tensor
 * \f[g_{ij}\equiv\braket{i}{j}\f]
 * is nondiagonal. The contravariant components of a state vector \f$\ket\Psi\f$ are then defined as the expansion coefficients
 * \f[\ket\Psi\equiv\Psi^i\ket{i},\f]
 * where we have adopted the convention that there is summation for indeces appearing twice. In this case one index has to be down, while the other one up, 
 * otherwise we have as it were a syntax error. This ensures that the result is independent of the choice of basis, as one would very well expect 
 * e.g. from the trace of an operator. The covariant components are the projections
 * \f[\Psi_i\equiv\braket{i}\Psi.\f]
 * The two are connected with the metric tensor: 
 * \f[\Psi_i=g_{ij}\Psi^j,\quad\Psi^i=g^{ij}\Psi_j,\f]
 * where it is easy to show that \f$g^{ij}\f$ is the matrix inverse of the above defined \f$g_{ij}\f$:
 * \f[g_{ik}g^{kj}\equiv g_i^j=\delta_{ij}.\f]
 * For the representation of operators assume \f$A\f$ is an operator, then
 * \f[\ket\Phi\equiv A\ket\Psi=\Psi^jA\ket{j},\f]
 * multiplying from the left by \f$\bra{i}\f$ yields
 * \f[\Phi_i=\bra iA\ket j\Psi^j\equiv A_{ij}\Psi^j.\f]
 * While the definition of the up-indeces matrix of the operator can be read from
 * \f[\Phi=\Phi^i\ket i\equiv A^{ij}\Psi_j\ket{i}=A^{ij}\braket{j}\Psi\ket{i}=\lp A^{ij}\ket{i}\bra{j}\rp\ket\Psi=A\ket\Psi.\f]
 * 
 * #### Hermitian conjugation
 * 
 * From the above it is easy to see that
 * \f[\lp A^\dag\rp_{ij}=\lp A_{ji}\rp^*,\;\text{and}\;\lp A^\dag\rp^{ij}=\lp A^{ji}\rp^*.\f]
 * The situation, however, is more involved with the mixed-index matrices since conjugation mixes the two kinds:
 * \f[{\lp A^\dag\rp_i}^j=\lp {A^j}_i\rp^*,\;\text{and}\;{\lp A^\dag\rp^i}_j=\lp {A_j}^i\rp^*.\f]
 * Here we prove the second:
 * \f[{A_j}^i=A_{jk}g^{ki}=\lp\lp A^\dag\rp_{kj}\rp^*g^{ki},\f]
 * so that
 * \f[\lp{A_j}^i\rp^*=\lp A^\dag\rp_{kj}\lp g^{ki}\rp^*=g^{ik}\lp A^\dag\rp_{kj}={\lp A^\dag\rp^i}_j,\f]
 * where we have used the Hermiticity of the metric tensor.
 * 
 * ### Implications
 * 
 * If one is to apply the above formalism for density matrices in nonorthogonal bases, one has to face some strange implications:
 * The well-known properties of Hermiticity and unit trace are split: \f$\rho_{ij}\f$ and \f$\rho^{ij}\f$ will be Hermitian matrices,
 * but their trace is irrelevant, in fact, it is not conserved during time evolution, 
 * while \f${\rho_i}^j\f$ and \f${\rho^i}_j\f$ have unit conserved trace.
 * 
 * This means that e.g. the quantumtrajectory::Master driver has to use either \f$\rho_{ij}\f$ or \f$\rho^{ij}\f$
 * (by convention, it uses the latter, which also means that quantumtrajectory::MCWF_Trajectory uses \f$\Psi^i\f$)
 * to be able to calculate the evolution using the real part as in \Eqref{eq:Master-MCWF} (otherwise, it is not the real part which appears there).
 * On the other hand, for all normalizations one index has to be pulled, and this (at least in the case of density operators)
 * has to be done in place to avoid excessive copying, which also means that the index has to be pulled back again afterwards.
 * 
 */


} // quantumdata

#endif // QUANTUMDATA_NONORTHOGONALSTATEVECTOR_H_INCLUDED
