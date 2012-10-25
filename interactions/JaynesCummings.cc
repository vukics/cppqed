#include "JaynesCummings_.h"

#include "impl/Pars.tcc"


namespace jaynescummings {


Pars::Pars(parameters::ParameterTable& p, const std::string& mod)
  : g(p.addMod("g",mod,"Qbit-mode coupling",dcomp(1.)))
{}


const structure::free::Tridiagonal sigmaop(qbit::Ptr qbit) {return qbit::sigmaop(qbit);}


const structure::free::Tridiagonal sigmaop(spin::Ptr spin) {return spin::sminus(spin);}


} // jaynescummings


