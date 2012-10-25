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

#include <boost/assign/list_of.hpp>


namespace jaynescummings {


struct Pars
{
  dcomp& g;

  operator dcomp&() const {return g;}

  Pars(parameters::ParameterTable&, const std::string& ="");

};

const structure::free::Tridiagonal sigmaop(qbit::Ptr);
const structure::free::Tridiagonal sigmaop(spin::Ptr);

class Base : public structure::Interaction<2>, public structure::TridiagonalHamiltonian<2,true>
{
protected:
  typedef structure::Interaction<2> IA_Base;
  typedef structure::TridiagonalHamiltonian<2,true> TDH_Base;

  template<typename QBIT_SPIN_BASE>
  Base(boost::shared_ptr<const QBIT_SPIN_BASE> qbitspin, mode::Ptr mode, const dcomp& g)
  : IA_Base(Frees(qbitspin,mode),RealFreqs(),boost::assign::tuple_list_of("g",g,sqrt(qbitspin->getDimension()*mode->getDimension()))),
    TDH_Base(tridiagMinusHC(conj(g)*jaynescummings::sigmaop(qbitspin)*mode::aop(mode).dagger())) {getParsStream()<<"# Jaynes-Cummings interaction\n";}

};


} // jaynescummings


#define BIG_NAMESPACE_NAME             jaynescummings
#define BIG_CLASS_NAME                 JaynesCummings
#define BIG_ADDITIONAL_PARAMETERS      , const dcomp& g
#define BIG_ADDITIONAL_PARAMETERS_PASS ,g

#include "details/BinaryInteractionGenerator.h"


#endif // ELEMENTS_INTERACTIONS_JAYNESCUMMINGS__H_INCLUDED
