// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_FREES_BICHROMATICMODE_H_INCLUDED
#define   CPPQEDELEMENTS_FREES_BICHROMATICMODE_H_INCLUDED

#include "Mode_.h"
#include "TridiagonalHamiltonian.h"


namespace mode {

struct ParsBichromatic : ParsPumpedLossy
{
  double &deltaOther;
  dcomp  &  etaOther;

  ParsBichromatic(parameters::Table&, const std::string& ="");
};

} // mode


template<bool TEMPERATURE=false, typename AveragingType=mode::Averaged>
class BichromaticMode
  : public mode::Liouvillean<TEMPERATURE>,
    public mode::Hamiltonian<true>,
    public ModeBase, public AveragingType
{
public:
  template<typename... AveragingConstructorParameters>
  BichromaticMode(const mode::ParsBichromatic& p, AveragingConstructorParameters&&... a)
    : mode::Liouvillean<TEMPERATURE>(p.kappa,p.nTh),
      mode::Hamiltonian<true>(0,dcomp(mode::finiteTemperatureHamiltonianDecay<TEMPERATURE>(p),-p.delta),p.eta,p.omegaKerr,p.omegaKerrAlter,p.cutoff),
      ModeBase(p.cutoff,RF{"deltaOther",p.deltaOther,1},{CF{"(kappa*(2*nTh+1),delta)",conj(get_zI()),1},CF{"eta",get_eta(),sqrt(p.cutoff)},CF{"etaOther",p.etaOther,sqrt(p.cutoff)}},"Bichromatic mode"),
      AveragingType(std::forward<AveragingConstructorParameters>(a)...),
      zI_Other(mode::finiteTemperatureHamiltonianDecay<TEMPERATURE>(p),-p.deltaOther)
  {
    mode::Hamiltonian<true>::getH_OverIs().push_back(
                                                     furnishWithFreqs(mode::pumping(p.etaOther,getDimension()),
                                                                      mode::mainDiagonal(zI_Other,0,getDimension()))
                                                     );
    getParsStream()<<"Bichromatic pumping."; mode::isFiniteTempStream<TEMPERATURE>(getParsStream(),p.nTh);
  }

private:
  const dcomp zI_Other;

};


namespace mode {


template<typename AveragingType, typename... AveragingConstructorParameters>
const Ptr make(const ParsBichromatic& p, QM_Picture qmp, AveragingConstructorParameters&&... a)
{
  if (!isNonZero(p.etaOther)) return make<AveragingType>(static_cast<const ParsPumpedLossy&>(p),qmp,a...);
  else {
    if (p.nTh) {return std::make_shared<BichromaticMode<true ,AveragingType> >(p,a...);}
    else       {return std::make_shared<BichromaticMode<false,AveragingType> >(p,a...);}
  }
}


} //mode


template class BichromaticMode<true>;
template class BichromaticMode<false>;


#endif // CPPQEDELEMENTS_FREES_BICHROMATICMODE_H_INCLUDED
