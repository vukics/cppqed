#include "JaynesCummings.h"

#include "Pars.h"

#include<boost/assign/list_of.hpp>

using namespace boost::assign;

using namespace mode;

namespace jaynescummings {


Pars::Pars(parameters::ParameterTable& p, const std::string& mod)
  : g(p.addMod("g",mod,"Qbit-mode coupling",dcomp(1.)))
{}


Base::Base(const QbitBase* qbit, const ModeBase* mode, const dcomp& g)
  : IA_Base(Frees(qbit,mode),RealFreqs(),tuple_list_of("g",g,sqrt(mode->getDimension()))),
    TDH_Base(tridiagMinusHC(conj(g)*qbit::sigmaop(qbit)*aop(mode).dagger()))
{
  getParsStream()<<"# Jaynes-Cummings interaction\n";
}

} // jaynescummings



/*
template<typename QbitType, typename ModeType>
JaynesCummings<QbitType,ModeType>::JaynesCummings(const QbitType& qbit, const ModeType& mode, const dcomp& g, mpl::bool_<IS_TD>)
  : JaynesCummingsBase(&qbit,&mode,g),
    TDH_Base(conj(g)*sigmaop()*aop(&mode).dagger()-g*sigmaop().dagger()*aop(&mode))
{
}

template<typename QbitType, typename ModeType>
JaynesCummings<QbitType,ModeType>::JaynesCummings(const QbitType& qbit, const ModeType& mode, const dcomp& g, mpl::bool_<IS_TD>)
  : JaynesCummingsBase(&qbit,&mode,g),
    TDH_Base(conj(g)*sigmaop()*aop(&mode).dagger()-g*sigmaop().dagger()*aop(&mode))
{
}
*/
