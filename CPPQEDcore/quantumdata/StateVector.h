// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#pragma once

#include "MultiArrayComplex.h"


namespace quantumdata {


using namespace ::cppqedutils;


using StorageType = multiarray::StorageType<dcomp>;


template<typename>
constexpr auto multiArrayRank_v=std::nullopt;

static const auto noInit = multiarray::noInit<dcomp>;

static const auto zeroInit = multiarray::zeroInit<dcomp>;


template <typename Derived>
struct VectorSpaceOperatorsGenerator
{

  friend auto operator+(const Derived& a, const Derived& b)
  {
    checkExtents(a,b,"VectorSpaceOperatorsGenerator plus");
    return Derived{getDimensions(a), [&] (size_t e) {
      auto r{noInit<multiArrayRank_v<Derived>>(e)}; for (auto& [re,ae,be] : std::views::zip(r,a.dataView,b.dataView)) re=ae+be; return r;
    }};
  }


  friend auto operator-(const Derived& a, const Derived& b)
  {
    checkExtents(a,b,"VectorSpaceOperatorsGenerator minus");
    return Derived{getDimensions(a), [&] (size_t e) {
      auto r{noInit<multiArrayRank_v<Derived>>(e)}; for (auto& [re,ae,be] : std::views::zip(r,a.dataView,b.dataView)) re=ae-be; return r;
    }};
  }


  friend auto operator*(dcomp v, const Derived& a)
  {
    return Derived{getDimensions(a), [&] (size_t e) {
      auto r{noInit<multiArrayRank_v<Derived>>(e)}; for (auto& [re,ae] : std::views::zip(r,a.dataView)) re=v*ae; return r;
    }};
  }

};


template <size_t RANK> struct DensityOperator; // forward declaration


template <size_t RANK>
using StateVectorView=MultiArrayView<dcomp,RANK>;


template <size_t RANK>
using StateVectorConstView=MultiArrayConstView<dcomp,RANK>;



/// State vector of arbitrary arity
/**
 * Owns its data. Should be present on the highest level of the simulations in a single copy
 * (representing the initial condition, to be evolved by the quantum dynamics).
 */
template<size_t RANK>
struct StateVector : MultiArray<dcomp,RANK>, private VectorSpaceOperatorsGenerator<StateVector<RANK>>
{
  using ABase = MultiArray<dcomp,RANK>;

  using Dimensions=Extents<RANK>;

  StateVector(const StateVector&) = delete; StateVector& operator=(const StateVector&) = delete;

  StateVector(StateVector&&) = default; StateVector& operator=(StateVector&&) = default;

  StateVector(MultiArray<dcomp,RANK>&& ma) : ABase{std::move(ma)} {}

  StateVector(Dimensions dimensions, auto&& initializer) : ABase{dimensions,std::forward<decltype(initializer)>(initializer)} {}

  explicit StateVector(Dimensions dimensions)
    : StateVector{dimensions, [] (size_t e) {auto res{zeroInit(e)}; res[0]=1.; return res;}} {}

  // StateVector(auto&& ... a) : ArrayBase<StateVector<RANK>>(std::forward<decltype(a)>(a)...) {}

  /// \name Metric
  //@{
    /// Returns the norm \f$\norm\Psi\f$, implemented in terms of ArrayBase::frobeniusNorm.
  friend double norm(const StateVector& sv) {return frobeniusNorm(sv);}
  
  friend double renorm(StateVector& sv)
  {
    double res=norm(sv);
    for (dcomp& v : sv.dataStorage()) v/=res;
    return res;
  }
  //@}
  
  /// Adds a dyad of the present object to `densityOperator`
  /**
   * This is done without actually forming the dyad in memory (so that this is not implemented in terms of StateVector::dyad),
   * which is important in situations when an average density operator is needed from an ensemble of state vectors, an example being quantumtrajectory::EnsembleMCWF. 
   */
  friend void addTo(const StateVector<RANK>& psi, DensityOperator<RANK>& rho, double weight=1.)
  {
    for (size_t size=psi.dataView.size(), i=0; i<size; ++i) for (size_t j=0; j<size; ++j)
      rho.dataStorage()[i+size*j]+=psi.dataView[i]*conj(psi.dataView[j]);
  }

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
  return directProduct(sv1, sv2, [] (dcomp v1, dcomp v2) {return v1*conj(v2);} );
}


/// Creates the direct product, relying on the direct-product constructor
template<size_t RANK1, size_t RANK2>
inline auto operator*(const StateVector<RANK1>& psi1, const StateVector<RANK2>& psi2)
{
  return StateVector<RANK1+RANK2>(directProduct<RANK1,RANK2>(psi1,psi2));
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
