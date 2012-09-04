// -*- C++ -*-
#ifndef ELEMENTS_INTERACTIONS_JAYNESCUMMINGS__H_INCLUDED
#define ELEMENTS_INTERACTIONS_JAYNESCUMMINGS__H_INCLUDED

// #include<boost/type_traits/is_base_of.hpp>

#include "JaynesCummingsFwd.h"

#include "Interaction.h"
#include "TridiagonalHamiltonian.h"

#include "Mode_.h"
#include "Qbit_.h"

#include "ParsFwd.h"


namespace jaynescummings {


struct Pars
{
  dcomp& g;

  Pars(parameters::ParameterTable&, const std::string& ="");

};


} // jaynescummings


class JaynesCummings : public structure::Interaction<2>, public structure::TridiagonalHamiltonian<2,true>
{
protected:
  typedef structure::Interaction<2> IA_Base;
  typedef structure::TridiagonalHamiltonian<2,true> TDH_Base;

  JaynesCummings(qbit::SmartPtr, mode::SmartPtr, const dcomp& g);

};





#endif // ELEMENTS_INTERACTIONS_JAYNESCUMMINGS__H_INCLUDED
