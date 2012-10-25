#include "JaynesCummings_.h"

#include "Qbit_.h"
#include "Spin.h"

#include "impl/Tridiagonal.tcc"

#include "impl/Pars.tcc"

#include <boost/assign/list_of.hpp>

using namespace boost::assign; using boost::dynamic_pointer_cast;

using namespace mode;

namespace jaynescummings {


Pars::Pars(parameters::ParameterTable& p, const std::string& mod)
  : g(p.addMod("g",mod,"Qbit-mode coupling",dcomp(1.)))
{}


Base::Base(QbitSpinCommonPtr spin, mode::Ptr mode, const dcomp& g)
  : IA_Base(Frees(spin,mode),RealFreqs(),tuple_list_of("g",g,sqrt(spin->getDimension()*mode->getDimension()))),
    TDH_Base(tridiagMinusHC(conj(g)*sigmaop(spin)*aop(mode).dagger()))
{
  getParsStream()<<"# Jaynes-Cummings interaction\n";
}


const structure::free::Tridiagonal sigmaop(QbitSpinCommonPtr qbitspin)
{
  if      (const qbit::Ptr qbit=dynamic_pointer_cast<const QbitBase>(qbitspin)) return sigmaop(qbit);
  else if (const spin::Ptr spin=dynamic_pointer_cast<const SpinBase>(qbitspin)) return spin::sminus (spin);
  else throw UnrecognizedPluginToJaynesCummings();
}

} // jaynescummings


