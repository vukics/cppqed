#include "JaynesCummings_.h"

#include "impl/Tridiagonal.tcc"

#include "impl/Pars.tcc"

#include <boost/assign/list_of.hpp>

using namespace boost::assign;

using namespace mode;


namespace jaynescummings {


Pars::Pars(parameters::ParameterTable& p, const std::string& mod)
  : g(p.addMod("g",mod,"Qbit-mode coupling",dcomp(1.)))
{}


} // jaynescummings


JaynesCummings::JaynesCummings(qbit::SmartPtr qbit, SmartPtr mode, const dcomp& g)
  : IA_Base(Frees(qbit,mode),RealFreqs(),tuple_list_of("g",g,sqrt(mode->getDimension()))),
    TDH_Base(tridiagMinusHC(conj(g)*qbit::sigmaop(qbit)*aop(mode).dagger()))
{
  getParsStream()<<"# Jaynes-Cummings interaction\n";
}




