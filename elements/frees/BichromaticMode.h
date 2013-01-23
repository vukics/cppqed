// -*- C++ -*-
#ifndef   ELEMENTS_FREES_BICHROMATICMODE_H_INCLUDED
#define   ELEMENTS_FREES_BICHROMATICMODE_H_INCLUDED

#include "BichromaticModeFwd.h"

#include "Mode_.h"
#include "impl/TridiagonalHamiltonian.tcc"

#include <boost/make_shared.hpp>


namespace mode {

struct ParsBichromatic : ParsPumpedLossy
{
  double &deltaOther;
  dcomp  &  etaOther;

  ParsBichromatic(parameters::ParameterTable&, const std::string& ="");
};


template<typename AveragingType, typename... AveragingConstructorParameters>
const Ptr make(const ParsBichromatic& p, QM_Picture qmp, const AveragingConstructorParameters&... a)
{
  using boost::make_shared;
  if (!isNonZero(p.etaOther)) return make<AveragingType>(static_cast<const ParsPumpedLossy&>(p),qmp,a...);
  else {
    if (p.nTh) {return make_shared<BichromaticMode<true ,AveragingType> >(p,a...);}
    else       {return make_shared<BichromaticMode<false,AveragingType> >(p,a...);}
  }
}


} //mode


template<bool IS_FINITE_TEMP, typename AveragingType>
class BichromaticMode
  : public mode::Liouvillean<IS_FINITE_TEMP>,
    public mode::Hamiltonian<true>,
    public ModeBase, public AveragingType
{
public:
  template<typename... AveragingConstructorParameters>
  BichromaticMode(const mode::ParsBichromatic& p, const AveragingConstructorParameters&... a)
    : mode::Liouvillean<IS_FINITE_TEMP>(p.kappa,p.nTh),
      mode::Hamiltonian<true>(0,dcomp(mode::finiteTemperatureHamiltonianDecay(p,*this),-p.delta),p.eta,p.cutoff),
      ModeBase(p.cutoff,
	       FREQS("deltaOther",p.deltaOther,1),
	       FREQS("(kappa*(2*nTh+1),delta)",conj(get_zI()),1)("eta",get_eta(),sqrt(p.cutoff))("etaOther",p.etaOther,sqrt(p.cutoff))
	       ),
      AveragingType(a...),
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


#endif // ELEMENTS_FREES_BICHROMATICMODE_H_INCLUDED
