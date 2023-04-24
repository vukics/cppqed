// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#pragma once

#include "MultiArrayComplex.h"


namespace quantumdata {

template<typename>
constexpr auto multiArrayRank_v=std::nullopt;

template <size_t RANK>
const auto noInit = ::cppqedutils::multiarray::noInit<dcomp,RANK>;

template <size_t RANK>
const auto zeroInit = ::cppqedutils::multiarray::zeroInit<dcomp,RANK>;

template <typename Derived>
struct VectorSpaceOperatorsGenerator
{

  friend Derived operator+(const Derived& a, const Derived& b)
  {
    checkExtents(a,b,"VectorSpaceOperatorsGenerator plus");
    Derived res{getDimensions(a), [&] (size_t e) {
      auto r{noInit<multiArrayRank_v<Derived>>(e)}; auto ri=r.begin();
      for (auto ai=a.dataView.begin(), bi=b.dataView.begin(); ai!=a.dataView.end(); *ri++=(*ai++)+(*bi++) ) ;
      return r;
    }};
    return res;
  }


  friend Derived operator-(const Derived& a, const Derived& b)
  {
    checkExtents(a,b,"VectorSpaceOperatorsGenerator minus");
    Derived res{getDimensions(a), [&] (size_t e) {
      auto r{noInit<multiArrayRank_v<Derived>>(e)}; auto ri=r.begin();
      for (auto ai=a.dataView.begin(), bi=b.dataView.begin(); ai!=a.dataView.end(); *ri++=(*ai++)-(*bi++) ) ;
      return r;
    }};
    return res;
  }


  friend Derived operator*(dcomp v, const Derived& a)
  {
    Derived res{getDimensions(a), [&] (size_t e) {
      auto r{noInit<multiArrayRank_v<Derived>>(e)}; auto ri=r.begin();
      for (auto ai=a.dataView.begin(); ai!=a.dataView.end(); (*ri++)=v*(*ai++) ) ;
      return r;
    }};
    return res;
  }

};


template <size_t RANK> struct DensityOperator; // forward declaration


template <size_t RANK>
using StateVectorView=cppqedutils::MultiArrayView<dcomp,RANK>;


template <size_t RANK>
using StateVectorConstView=cppqedutils::MultiArrayConstView<dcomp,RANK>;


/// State vector of arbitrary arity
/**
 * Owns its data. Should be present on the highest level of the simulations in a single copy
 * (representing the initial condition, to be evolved by the quantum dynamics).
 */
template<size_t RANK>
struct StateVector : ::cppqedutils::MultiArray<dcomp,RANK>, private VectorSpaceOperatorsGenerator<StateVector<RANK>>
{
  using ABase = ::cppqedutils::MultiArray<dcomp,RANK>;

  using Dimensions=::cppqedutils::Extents<RANK>;

  StateVector(const StateVector&) = delete; StateVector& operator=(const StateVector&) = delete;

  StateVector(StateVector&&) = default; StateVector& operator=(StateVector&&) = default;

  StateVector(::cppqedutils::MultiArray<dcomp,RANK>&& ma) : ABase{std::move(ma)} {}

  StateVector(Dimensions dimensions, auto&& initializer) : ABase{dimensions,std::forward<decltype(initializer)>(initializer)} {}

  explicit StateVector(Dimensions dimensions)
    : StateVector{dimensions, [] (size_t e) {auto res{zeroInit<RANK>(e)}; res[0]=1.; return res;}} {}

  // StateVector(auto&& ... a) : ArrayBase<StateVector<RANK>>(std::forward<decltype(a)>(a)...) {}

  /// \name Metric
  //@{
    /// Returns the norm \f$\norm\Psi\f$, implemented in terms of ArrayBase::frobeniusNorm.
  friend double norm(const StateVector& sv) {return frobeniusNorm(sv);}
  
  friend double renorm(StateVector& sv) {double res=norm(sv); std::ranges::for_each(sv.mutableView().dataView,[=] (auto& v) {v/=res;}); return res;}
  //@}
  
  /// Adds a dyad of the present object to `densityOperator`
  /**
   * This is done without actually forming the dyad in memory (so that this is not implemented in terms of StateVector::dyad),
   * which is important in situations when an average density operator is needed from an ensemble of state vectors, an example being quantumtrajectory::EnsembleMCWF. 
   */
  void addTo(DensityOperator<RANK>& rho, double weight=1.) const;/*
  {
    using namespace linalg;
    CMatrix matrix(rho.matrixView());
    CVector vector(vectorView());
    size_t dim(this->getTotalDimension());
    for (size_t i=0; i<dim; i++) for (size_t j=0; j<dim; j++) matrix(i,j)+=weight*vector(i)*conj(vector(j));
  }*/
  /*
  friend Derived operator+(const Derived&, const Derived&);
  friend Derived operator-(const Derived& aa, const Derived& bb)
  {
    auto res{Derived::clone(aa)};
    for (auto&& [r,a,b] : boost::combine(res.mutableView().dataView,aa.dataView,bb.dataView) ) r=a-b;
    return res;
  };

  friend Derived operator*(dcomp, const ArrayBase&);
*/
};


template <size_t RANK>
auto getDimensions(const StateVector<RANK>& psi) {return psi.extents;}


/// Forms a dyad with the argument, a rather expensive operation, implemented as a directProduct
template <size_t RANK>
auto dyad(const StateVector<RANK>& sv1, const StateVector<RANK>& sv2)
{
#ifndef NDEBUG
  if (sv1.extents != sv2.extents) throw std::runtime_error("Mismatch in StateVector::dyad dimensions");
#endif // NDEBUG
  return cppqedutils::directProduct(sv1, sv2, [] (dcomp v1, dcomp v2) {return v1*conj(v2);} );
}


/// Creates the direct product, relying on the direct-product constructor
template<size_t RANK1, size_t RANK2>
inline auto operator*(const StateVector<RANK1>& psi1, const StateVector<RANK2>& psi2)
{
  return StateVector<RANK1+RANK2>(cppqedutils::directProduct<RANK1,RANK2>(psi1,psi2));
}


/// Calculates the inner product, relying on StateVector::vectorView
template<size_t RANK>
dcomp braket(const StateVector<RANK>& psi1, const StateVector<RANK>& psi2);/*
{
  using blitz::tensor::i;
  linalg::CVector temp(psi1.getTotalDimension());
  temp=conj(psi1.vectorView()(i))*psi2.vectorView()(i);
  return sum(temp);
}*/


template <size_t RANK>
constexpr auto multiArrayRank_v<StateVector<RANK>> = RANK;


} // quantumdata


template <size_t RANK>
constexpr auto cppqedutils::passByValue_v<quantumdata::StateVector<RANK>> = false;


namespace structure {

using ::quantumdata::StateVectorView, ::quantumdata::StateVectorConstView;

} // structure


/*
namespace boost { namespace numeric { namespace odeint {

template<typename Derived>
struct is_resizeable<::quantumdata::ArrayBase<Derived>> : boost::true_type {};

template<typename Derived>
struct same_size_impl<::quantumdata::ArrayBase<Derived>, ::quantumdata::ArrayBase<Derived>>
{ // define how to check size
  static bool same_size(const ::quantumdata::ArrayBase<Derived> &v1, const ::quantumdata::ArrayBase<Derived> &v2) {return v1.extents == v2.extents;}
};

/// TODO: a reserve could be defined for the vector to be resized
template<typename Derived>
struct resize_impl<::quantumdata::ArrayBase<Derived>, ::quantumdata::ArrayBase<Derived>>
{ // define how to resize
  static void resize(::quantumdata::ArrayBase<Derived> &v1, const ::quantumdata::ArrayBase<Derived> &v2) {v1.resize( v2.extents );}
};

template<typename Derived>
struct vector_space_norm_inf<::quantumdata::ArrayBase<Derived>>
{
  typedef double result_type;
  double operator()(const ::quantumdata::ArrayBase<Derived>& v ) const
  {
    return max( abs(v) );
  }
};

template<typename Derived>
struct norm_result_type<::quantumdata::ArrayBase<Derived>> : mpl::identity<double> {};

} } } // boost::numeric::odeint
*/
