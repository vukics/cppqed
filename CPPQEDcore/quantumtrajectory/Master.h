// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "QuantumSystemDynamics.h"
#include "QuantumTrajectory.h"



namespace quantumtrajectory {


/// Auxiliary tools to Master
namespace master {


typedef cppqedutils::ode::Pars<> Pars;


} // master


/// Master equation evolution from a \link quantumdata::DensityOperator density-operator\endlink initial condition
/**
 * \see \ref masterequation.
 * 
 * \note The ODE driver underlying this class needs to store several (typically 6–7, this depends on the chosen driver) density-operator instants.
 *
 */
template <size_t RANK,
          ::structure::quantum_system_dynamics<RANK> QSD,
          ::cppqedutils::ode::engine<typename ::quantumdata::DensityOperator<RANK>::StorageType> OE>
struct Master
{
  typedef quantumdata::DensityOperator<RANK> DensityOperator;

  Master(auto&& qsd, auto&& rho, auto&& oe)
    : qsd{std::forward<decltype(qsd)>(qsd)}, rho{std::forward<decltype(rho)>(rho)}, oe{std::forward<decltype(oe)>(oe)},
      rowIterationOffsets_{::cppqedutils::calculateSlicesOffsets<rowIterationRetainedAxes>(rho.extents)} {}

  double time=0., time0=0.;
  QSD qsd;
  DensityOperator rho;
  OE oe;

  friend double getDtDid(const Master& m) {return getDtDid(m.oe);}

  friend double getTime(const Master& m) {return m.time;}

  /// TODO: put here the system-specific things
  friend ::cppqedutils::LogTree logIntro(const Master& m)
  {
    return {{"Master",{{"odeEngine",logIntro(m.oe)},{"System","TAKE FROM SYSTEM"}}}};
    // streamCharacteristics(m.qsd,streamParameters(m.qsd,streamIntro(m.oe,os)<<"Solving Master equation."<<endl<<endl ) )<<endl;
  }

  friend ::cppqedutils::LogTree logOutro(const Master& m) {return logOutro(m.oe);}

  friend ::cppqedutils::LogTree step(Master& m, double deltaT)
  {
    // TODO: think over what happens here when there is no Hamiltonian to apply – it’s probably OK, since at least one superoperator is always there
    // Moreover, I’m not sure anymore that these algorithms should be prepared for such special cases.
    auto res = step ( m.oe, deltaT, [&] (const typename DensityOperator::StorageType& rhoRaw,
                                         typename DensityOperator::StorageType& drhodtRaw, double t)
    {
      ::quantumdata::DensityOperatorConstView<RANK> rho{m.rho.extents,m.rho.strides,0,rhoRaw};
      ::quantumdata::DensityOperatorView<RANK> drhodt{m.rho.extents,m.rho.strides,0,drhodtRaw.begin(),std::ranges::fill(drhodtRaw,0)};

      // TODO: std::ranges::views::zip to be applied here, but for some reason, sliceRange is not compatible with zipping – it does work with sliceRangeSimple, though.
      //for (auto [psi,dpsidt] : std::views::zip(::cppqedutils::sliceRange<Master::rowIterationRetainedAxes>(rho,m.rowIterationOffsets_),
      //                                         ::cppqedutils::sliceRange<Master::rowIterationRetainedAxes>(drhodt,m.rowIterationOffsets_) ) )
      //  ::structure::applyHamiltonian(getHa(m.qsd),t,psi,dpsidt,m.time0) ;
      {
        auto psiRange{::cppqedutils::sliceRange<Master::rowIterationRetainedAxes>(rho,m.rowIterationOffsets_)};
        auto dpsidtRange{::cppqedutils::sliceRange<Master::rowIterationRetainedAxes>(drhodt,m.rowIterationOffsets_)};
        for (auto [psi,dpsidt]=std::make_tuple(psiRange.begin(),dpsidtRange.begin()); psi!=psiRange.end();
             ::structure::applyHamiltonian(getHa(m.qsd),t,*psi++,*dpsidt++,m.time0)) ;
      }

      twoTimesRealPartOfSelf(drhodt);

      for (auto&& lindblad : getLi(m.qsd))
        ::structure::applySuperoperator<RANK>(lindblad.superoperator, t, rho, drhodt);

    },m.time,m.rho.dataStorage());

    // exact propagation
    for (auto&& psi : ::cppqedutils::sliceRange<Master::rowIterationRetainedAxes>(m.rho.mutableView(),m.rowIterationOffsets_))
      ::structure::applyPropagator(getHa(m.qsd),m.time,psi,m.time0);
    hermitianConjugateSelf(m.rho);
    for (auto&& psi : ::cppqedutils::sliceRange<Master::rowIterationRetainedAxes>(m.rho.mutableView(),m.rowIterationOffsets_))
      ::structure::applyPropagator(getHa(m.qsd),m.time,psi,m.time0);
    conj(m.rho);
    m.time0=m.time;

    // The following "smoothening" of rho_ has proven necessary for the algorithm to remain stable:
    // We make the approximately Hermitian and normalized rho_ exactly so.

    twoTimesRealPartOfSelf(m.rho.mutableView());

    renorm(m.rho);

    return res;
  }

  friend ::cppqedutils::LogTree dataStreamKey(const Master& m) {return {{"Master","TAKE FROM SYSTEM"}};}

  friend auto temporalDataPoint(const Master& m)
  {
    return ::structure::calculateAndPostprocess<RANK>( getEV(m.qsd), m.time, ::quantumdata::DensityOperatorConstView<RANK>(m.rho) );
  }

  friend ::cppqedutils::iarchive& readFromArrayOnlyArchive(Master& m, ::cppqedutils::iarchive& iar) {return iar & m.rho;} // MultiArray can be (de)serialized

  /** structure of Master archives:
  * metaData – array – time – ( odeStepper – odeLogger – dtDid – dtTry )
  * state should precede time in order to be compatible with array-only archives
  */
  template <typename Archive>
  friend auto& stateIO(Master& m, Archive& ar)
  {
    stateIO(m.oe, ar & m.rho & m.time);
    // The following is a no-op in the case of state output, since time0 is equated with time @ the end of each step, however, it is an essential step for state input!
    m.time0=m.time;
    return ar;
  }


private:
  const std::vector<size_t> rowIterationOffsets_;

  static constexpr auto rowIterationRetainedAxes=::cppqedutils::compileTimeOrdinals<RANK>;

  // const DensityOperatorStreamer<RANK,V> dos_;

};


template<typename QSD, typename DO, typename OE>
Master(QSD , DO , OE ) -> Master< quantumdata::multiArrayRank_v<DO>/2, QSD, OE >;


/*
namespace master {

template<typename ODE_Engine, typename V, typename StateVector_OR_DensityOperator>
auto make(::structure::QuantumSystemPtr<std::decay_t<StateVector_OR_DensityOperator>::N_RANK> sys,
          StateVector_OR_DensityOperator&& state, const Pars& p, EntanglementMeasuresSwitch ems)
{
  return Master<std::decay_t<StateVector_OR_DensityOperator>::N_RANK,ODE_Engine,V>(
    sys,std::forward<StateVector_OR_DensityOperator>(state),ODE_Engine{initialTimeStep(sys),p},ems);
}
  
} // master

*/
} // quantumtrajectory


template <size_t RANK, typename QSD, typename OE>
struct cppqedutils::trajectory::MakeSerializationMetadata<quantumtrajectory::Master<RANK,QSD,OE>>
{
  static auto _() {return SerializationMetadata{"CArray","Master",2*RANK};}
};


/*


template<int RANK, typename ODE_Engine, typename V>
std::ostream& quantumtrajectory::Master<RANK,ODE_Engine,V>::streamParameters(std::ostream& os) const
{
  using namespace std;

  ::structure::streamCharacteristics(sys_, sys_->streamParameters(ode_.streamParameters(os)<<"Solving Master equation."<<endl<<endl ) )<<endl;

  if (const auto li=::structure::castLi(sys_)) {
    os<<"Decay channels:\n";
    {
      size_t i=0;
      li->streamKey(os,i);
    }
    os<<"Explicit superoperator calculations: ";
    DensityOperator rhotemp(rho_.getDimensions());
    {
      int n=0;
      for (int i=0; i<li->nAvr(); ++i)
        try {
          li->actWithSuperoperator(0,rho_.getArray(),rhotemp.getArray(),i);
          os<<i<<' '; ++n;
        }
        catch (const ::structure::SuperoperatorNotImplementedException&) {}
      if (!n) os<<"none";
    }
    os<<endl;

  }
  
  return os;
}

*/

/*

namespace boost { namespace numeric { namespace odeint {

template <size_t RANK>
struct is_resizeable<::quantumdata::DensityOperator<RANK>> : boost::true_type {};

template <size_t RANK>
struct same_size_impl<::quantumdata::DensityOperator<RANK>, ::quantumdata::DensityOperator<RANK>>
{ // define how to check size
  static bool same_size(const ::quantumdata::DensityOperator<RANK> &v1, const ::quantumdata::DensityOperator<RANK> &v2) {return v1.extents == v2.extents;}
};

/// TODO: a reserve could be defined for the vector to be resized
template <size_t RANK>
struct resize_impl<::quantumdata::DensityOperator<RANK>, ::quantumdata::DensityOperator<RANK>>
{ // define how to resize
  static void resize(::quantumdata::DensityOperator<RANK> &v1, const ::quantumdata::DensityOperator<RANK> &v2) {v1.resize( v2.extents );}
};

template <size_t RANK>
struct vector_space_norm_inf<::quantumdata::DensityOperator<RANK>>
{
  typedef double result_type;
  double operator()(const ::quantumdata::DensityOperator<RANK>& v ) const
  {
    return max( abs(v) );
  }
};

template <size_t RANK>
struct norm_result_type<::quantumdata::DensityOperator<RANK>> : mpl::identity<double> {};

} } } // boost::numeric::odeint

*/
