// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "JaynesCummings_.h"

#include "Pars.h"


namespace jaynescummings {


Pars::Pars(parameters::Table& p, const std::string& mod)
  : g(p.add("g",mod,"Qbit-mode coupling",dcomp(1.)))
{}


const structure::freesystem::Tridiagonal sigmaop(qbit::Ptr qbit) {return qbit::sigmaop(qbit);}


const structure::freesystem::Tridiagonal sigmaop(spin::Ptr spin) {return spin::sminus(spin);}


} // jaynescummings


