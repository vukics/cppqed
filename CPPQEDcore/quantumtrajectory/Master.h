// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

// #include "QuantumTrajectory.h"

#include "QuantumSystemDynamics.h"

#include "Trajectory.h"



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
          ::cppqedutils::ode::engine<::quantumdata::DensityOperator<RANK>> OE>
struct Master
{
  typedef quantumdata::DensityOperator<RANK> DensityOperator;

  friend double getDtDid(const Master& m) {return getDtDid(m.oe);}

  friend double getTime(const Master& m) {return m.time;}

  /// TODO: put here the system-specific things
  friend ::cppqedutils::LogTree logIntro(const Master& m)
  {
    return {{"Master",{{"odeEngine",logIntro(m.oe)},{"System","TAKE FROM SYSTEM"}}}};
    // streamCharacteristics(m.qsd,streamParameters(m.qsd,streamIntro(m.oe,os)<<"Solving Master equation."<<endl<<endl ) )<<endl;
  }

  friend ::cppqedutils::LogTree logOutro(const Master& m) {return logOutro(m.oe);}

  friend ::cppqedutils::LogTree step(Master& m, double deltaT) {
    auto res = step ( m.oe, deltaT, [&] (const ::quantumdata::DensityOperator<RANK>& rho, ::quantumdata::DensityOperator<RANK>& drhodt, double t)
    {
      for (dcomp& v : drhodt.dataView) v=0;

      for (auto&& [psi,dpsidt] : boost::combine(::cppqedutils::sliceRange<Master::rowIterationRetainedAxes>(rho,m.rowIterationOffsets),
                                                ::cppqedutils::sliceRange<Master::rowIterationRetainedAxes>(drhodt.mutableView(),m.rowIterationOffsets)))
        getHa(m.qsd)(t,psi,dpsidt,m.time0);

      twoTimesRealPartOfSelf(drhodt);

      for (auto&& lindblad : getLi(m.qsd))
        applySuperoperator(lindblad.superoperator, t, rho, drhodt.mutableView());

    },m.time,m.rho.mutableView());

    // exact propagation
    for (auto&& psi : ::cppqedutils::sliceRange<Master::rowIterationRetainedAxes>(m.rho,m.rowIterationOffsets)) getHa(m.qsd)(m.time,psi,m.time0);
    hermitianConjugateSelf(m.rho);
    for (auto&& psi : ::cppqedutils::sliceRange<Master::rowIterationRetainedAxes>(m.rho,m.rowIterationOffsets)) getHa(m.qsd)(m.time,psi,m.time0);
    conj(m.rho);
    m.time0=m.time;

    // The following "smoothening" of rho_ has proven necessary for the algorithm to remain stable:
    // We make the approximately Hermitian and normalized rho_ exactly so.

    twoTimesRealPartOfSelf(m.rho);

    renorm(m.rho);

    return res;
  }

  friend ::cppqedutils::LogTree dataStreamKey(const Master& m) {return {{"Master","TAKE FROM SYSTEM"}};}

  friend auto temporalDataPoint(const Master& m) { return calculateAndPostprocess( getEv(m.qsd), m.time, m.rho );  }

  friend ::cppqedutils::iarchive& readFromArrayOnlyArchive(Master& m, ::cppqedutils::iarchive& iar) {return iar & m.rho;} // MultiArray can be (de)serialized

  /** structure of Master archives:
  * metaData – array – time – ( odeStepper – odeLogger – dtDid – dtTry )
  * state should precede time in order to be compatible with array-only archives
  */
  template <typename Archive>
  friend auto& stateIO(Master& m, Archive& ar)
  {
    stateIO(m.ode, ar & m.rho & m.time);
    // The following is a no-op in the case of state output, since time0 is equated with time @ the end of each step, however, it is an essential step for state input!
    m.time0=m.time;
    return ar;
  }


  double time=0., time0=0.;
  QSD qsd;
  DensityOperator rho;
  OE oe;

private:
  const std::vector<size_t> rowIterationOffsets;

  static constexpr auto rowIterationRetainedAxes=hana::range_c<size_t,0,RANK>;

  // const DensityOperatorStreamer<RANK,V> dos_;

};


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

/*
template <int RANK, typename ODE_Engine, typename V>
struct cppqedutils::trajectory::MakeSerializationMetadata<quantumtrajectory::Master<RANK,ODE_Engine,V>>
{
  static auto _() {return SerializationMetadata{"CArray","Master",2*RANK};}
};




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
