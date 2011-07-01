#include "Qbit.h"

#include "ParsQbit.h"

#include<boost/assign/list_of.hpp>

#include<boost/bind.hpp>


using namespace std;
using namespace boost;
using namespace assign;
using namespace mathutils;


namespace qbit {


///////////
//
// Averaged
//
///////////

Averaged::Averaged() 
  : Base("Qbit",list_of("rho00")("rho11")("real(rho10=<sigma>)")("imag(\")"))
{
}

const Averaged::Averages Averaged::average(const LazyDensityOperator& matrix) const
{
  Averages averages(4);
  averages(0)=matrix(0);
  averages(1)=matrix(1);
  averages(2)=real(matrix(1,0));
  averages(3)=imag(matrix(1,0));
  return averages;
}


} // qbit

////////////////
//
// Highest level
//
////////////////


QbitBase::QbitBase(const RealFreqs& realFreqs, const ComplexFreqs& complexFreqs)
  : ModeBase(2,realFreqs,complexFreqs), Averaged()
{
  getParsStream()<<"# Qbit\n";
}

#define BASE_init(R,C) QbitBase((R),(C))
#define BASE_initR(R)  BASE_init((R),ComplexFreqs())
#define BASE_initC(C)  BASE_init(RealFreqs(),(C))

#define TUPLE_delta(ISIP) "delta",-imag(BOOST_PP_IF(ISIP,get_zI,get_zSch)()),1//,BOOST_PP_IF(ISIP,1,p.cutoff)
#define TUPLE_eta "eta",get_eta(),1//,sqrt(p.cutoff)
#define TUPLE_gammadelta(ISIP) "(gamma,delta)",conj(BOOST_PP_IF(ISIP,get_zI,get_zSch)()),1//,BOOST_PP_IF(ISIP,1,p.cutoff)
#define TUPLE_gamma "gamma",real(get_zSch()),1//,p.cutoff

Qbit::Qbit(const qbit::Pars& p)
  : Exact(dcomp(0,-p.delta)),
    BASE_initR(tuple_list_of(TUPLE_delta(1)))
{}


QbitSch::QbitSch(const qbit::Pars& p)
  : qbit::Hamiltonian<false>(dcomp(0,-p.delta),0),
    BASE_initR(tuple_list_of(TUPLE_delta(0)))
{
  getParsStream()<<"# Schroedinger picture.\n";
}

   
PumpedQbit::PumpedQbit(const qbit::ParsPumped& p)
  : qbit::Hamiltonian<true>(0,dcomp(0,-p.delta),p.eta),
    BASE_init(tuple_list_of(TUPLE_delta(1)),tuple_list_of(TUPLE_eta))
{
  getParsStream()<<"# Pumped.\n";
}


PumpedQbitSch::PumpedQbitSch(const qbit::ParsPumped& p)
  : qbit::Hamiltonian<false>(dcomp(0,-p.delta),p.eta),
    BASE_init(tuple_list_of(TUPLE_delta(0)),tuple_list_of(TUPLE_eta))
{
  getParsStream()<<"# Pumped, Schroedinger picture.\n";
}


LossyQbit::LossyQbit(const qbit::ParsLossy& p)
  : qbit::Liouvillean(p.gamma),
    Exact(dcomp(p.gamma,-p.delta)),
    BASE_initC(tuple_list_of(TUPLE_gammadelta(1)))
{
  getParsStream()<<"# Lossy.\n";
}


LossyQbitSch::LossyQbitSch(const qbit::ParsLossy& p)
  : qbit::Liouvillean(p.gamma),
    qbit::Hamiltonian<false>(dcomp(p.gamma,-p.delta),0),
    BASE_initR(tuple_list_of(TUPLE_gamma)(TUPLE_delta(0)))
{
  getParsStream()<<"# Lossy, Schroedinger picture.\n";
}


LossyQbitUIP::LossyQbitUIP(const qbit::ParsLossy& p)
  : qbit::Liouvillean(p.gamma),
    qbit::Hamiltonian<true>(dcomp(p.gamma,0),dcomp(0,-p.delta),0),
    BASE_initR(tuple_list_of(TUPLE_gamma)(TUPLE_delta(1)))
{
  getParsStream()<<"# Lossy, Unitary interaction picture.\n";
}


PumpedLossyQbit::PumpedLossyQbit(const qbit::ParsPumpedLossy& p)
  : qbit::Liouvillean(p.gamma),
    qbit::Hamiltonian<true>(0,dcomp(p.gamma,-p.delta),p.eta),
    BASE_initC(tuple_list_of(TUPLE_gammadelta(1))(TUPLE_eta))
{
  getParsStream()<<"# Lossy, Unitary interaction picture.\n";
}


PumpedLossyQbitUIP::PumpedLossyQbitUIP(const qbit::ParsPumpedLossy& p)
  : qbit::Liouvillean(p.gamma),
    qbit::Hamiltonian<true>(dcomp(p.gamma,0),dcomp(0,-p.delta),p.eta),
    BASE_init(tuple_list_of(TUPLE_gamma)(TUPLE_delta(1)),tuple_list_of(TUPLE_eta))
{
  getParsStream()<<"# PumpedLossy, Unitary interaction picture.\n";
}


PumpedLossyQbitSch::PumpedLossyQbitSch(const qbit::ParsPumpedLossy& p)
  : qbit::Liouvillean(p.gamma),
    qbit::Hamiltonian<false>(dcomp(p.gamma,-p.delta),p.eta),
    BASE_initC(tuple_list_of(TUPLE_gammadelta(0))(TUPLE_eta))
{
  getParsStream()<<"# PumpedLossy, Schroedinger picture.\n";
}

#undef  TUPLE_gamma
#undef  TUPLE_gammadelta
#undef  TUPLE_eta
#undef  TUPLE_delta

#undef  BASE_initR
#undef  BASE_initC
#undef  BASE_init

//////////
//
// Helpers
//
//////////

namespace qbit {


SmartPtr maker(const ParsPumpedLossy& p, QM_Picture qmp)
{
  switch (qmp) {
  case QMP_IP  :
    return SmartPtr(new PumpedLossyQbit   (p));
  case QMP_UIP :
    return SmartPtr(new PumpedLossyQbitUIP(p));
  case QMP_SCH :
    ;
  }
  return SmartPtr(new PumpedLossyQbitSch(p));
}


} // qbit
