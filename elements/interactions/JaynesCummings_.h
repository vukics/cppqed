// -*- C++ -*-
#ifndef ELEMENTS_INTERACTIONS_JAYNESCUMMINGS__H_INCLUDED
#define ELEMENTS_INTERACTIONS_JAYNESCUMMINGS__H_INCLUDED

// #include<boost/type_traits/is_base_of.hpp>

#include "JaynesCummingsFwd.h"

#include "Interaction.h"
#include "TridiagonalHamiltonian.h"
#include "impl/Tridiagonal.tcc"

#include "Mode_.h"
#include "Qbit_.h"
#include "Spin.h"

#include "ParsFwd.h"

#include <boost/make_shared.hpp>


namespace jaynescummings {

typedef boost::shared_ptr<const Base<false> > Ptr;

struct Pars
{
  dcomp& g;

  operator dcomp&() const {return g;}

  Pars(parameters::ParameterTable&, const std::string& ="");

};

const structure::free::Tridiagonal sigmaop(qbit::Ptr);
const structure::free::Tridiagonal sigmaop(spin::Ptr);


template<>
class Base<false> : public structure::Interaction<2>
{
protected:
  typedef structure::Interaction<2> IA_Base;

  template<typename QBIT_SPIN_BASE>
  Base(boost::shared_ptr<const QBIT_SPIN_BASE> qbitspin, mode::Ptr mode, const dcomp& g=0)
  : IA_Base(Frees(qbitspin,mode),RealFreqs(),FREQS("g",g,sqrt(qbitspin->getDimension()*mode->getDimension())))
  {
    getParsStream()<<"# Jaynes-Cummings interaction\n";
  }

};


template<>
class Base<true> : public Base<false>, public structure::TridiagonalHamiltonian<2,true>
{
protected:
  typedef structure::TridiagonalHamiltonian<2,true> TDH_Base;

  template<typename QBIT_SPIN_BASE>
  Base(boost::shared_ptr<const QBIT_SPIN_BASE> qbitspin, mode::Ptr mode, const dcomp& g)
  : Base<false>(qbitspin,mode,g),
    TDH_Base(tridiagMinusHC(conj(g)*jaynescummings::sigmaop(qbitspin)*mode::aop(mode).dagger()))
  {}

};


} // jaynescummings


#define BIG_NAMESPACE_NAME                      jaynescummings
#define BIG_CLASS_NAME                          JaynesCummings
#define BIG_ADDITIONAL_PARAMETERS               , const dcomp& g
#define BIG_ADDITIONAL_PARAMETERS_PASS          ,g
#define BIG_ADDITIONAL_TEMPLATE_PARAMETERS      bool IS_HA,
#define BIG_ADDITIONAL_TEMPLATE_PARAMETERS_PASS <IS_HA>

#include "details/BinaryInteractionGenerator.h"


namespace jaynescummings {


template<typename A, typename F1, typename F2>
const Ptr make(const F1& f1, const F2& f2, const dcomp& g, const A& =A())
{
  if (isNonZero(g)) return boost::make_shared<JaynesCummings<true ,A> >(f1,f2,g );
  else              return boost::make_shared<JaynesCummings<false,A> >(f1,f2,0.);
}


template<typename F1, typename F2>
const Ptr make(const F1& f1, const F2& f2, const dcomp& g)
{
  return make<EmptyAveragingBaseForInteractions>(f1,f2,g);  
}


} // jaynescummings


#endif // ELEMENTS_INTERACTIONS_JAYNESCUMMINGS__H_INCLUDED
