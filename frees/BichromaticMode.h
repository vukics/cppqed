// -*- C++ -*-
#ifndef   BICHROMATIC_MODE_ELEMENT_INCLUDED
#define   BICHROMATIC_MODE_ELEMENT_INCLUDED

#include "BichromaticModeFwd.h"

#include "Mode.h"


namespace mode {

struct ParsBichromatic : ParsPumpedLossy
{
  double &deltaOther;
  dcomp  &  etaOther;

  ParsBichromatic(parameters::ParameterTable&, const std::string& ="");
};


template<typename A>
const SmartPtr maker(const ParsBichromatic& p, QM_Picture qmp, const A& a)
{
  if (!isNonZero(p.etaOther)) return maker(static_cast<const ParsPumpedLossy&>(p),qmp,a);
  else {
    if (p.nTh) {return SmartPtr(new BichromaticMode<true ,A>(p,a));}
    else       {return SmartPtr(new BichromaticMode<false,A>(p,a));}
  }
}


} //mode


template<bool IS_FINITE_TEMP, typename A>
class BichromaticMode
  : public mode::Liouvillean<IS_FINITE_TEMP>,
    public mode::Hamiltonian<true>,
    public ModeBase, public A
{
public:
  BichromaticMode(const mode::ParsBichromatic& p, const A& a)
    : mode::Liouvillean<IS_FINITE_TEMP>(p.kappa,p.nTh),
      mode::Hamiltonian<true>(0,dcomp(mode::finiteTemperatureHamiltonianDecay(p,*this),-p.delta),p.eta,p.cutoff),
      ModeBase(p.cutoff,
	       boost::assign::tuple_list_of("deltaOther",p.deltaOther,1),
	       boost::assign::tuple_list_of("(kappa*(2*nTh+1),delta)",conj(get_zI()),1)("eta",get_eta(),sqrt(p.cutoff))("etaOther",p.etaOther,sqrt(p.cutoff))
	       ),
      A(a),
      zI_Other(mode::finiteTemperatureHamiltonianDecay(p,*this),-p.deltaOther)
  {
    mode::Hamiltonian<true>::getH_OverIs().push_back(
						     furnishWithFreqs(mode::pumping(p.etaOther,getDimension()),
								      mode::mainDiagonal(zI_Other,getDimension()))
						     );
    getParsStream()<<"# Bichromatic pumping."; mode::isFiniteTempStream(getParsStream(),p.nTh,boost::mpl::bool_<IS_FINITE_TEMP>());
  }

private:
  const dcomp zI_Other;

};



template class BichromaticMode<true>;
template class BichromaticMode<false>;


#endif // BICHROMATIC_MODE_ELEMENT_INCLUDED
