// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDELEMENTS_FREES_BICHROMATICMODE_H_INCLUDED
#define   CPPQEDELEMENTS_FREES_BICHROMATICMODE_H_INCLUDED

#include "BichromaticModeFwd.h"

#include "Mode_.h"
#include "TridiagonalHamiltonian.tcc"

#include <boost/make_shared.hpp>


namespace mode {

struct ParsBichromatic : ParsPumpedLossy
{
  double &deltaOther;
  dcomp  &  etaOther;

  ParsBichromatic(parameters::ParameterTable&, const std::string& ="");
};


template<typename AveragingType, typename... AveragingConstructorParameters>
const Ptr make(const ParsBichromatic& p, QM_Picture qmp, AveragingConstructorParameters&&... a)
{
  if (!isNonZero(p.etaOther)) return make<AveragingType>(static_cast<const ParsPumpedLossy&>(p),qmp,a...);
  else {
    if (p.nTh) {return boost::make_shared<BichromaticMode<true ,AveragingType> >(p,a...);}
    else       {return boost::make_shared<BichromaticMode<false,AveragingType> >(p,a...);}
  }
}


} //mode


template<bool TEMPERATURE, typename AveragingType>
class BichromaticMode
  : public mode::Liouvillean<TEMPERATURE>,
    public mode::Hamiltonian<true>,
    public ModeBase, public AveragingType
{
public:
  template<typename... AveragingConstructorParameters>
  BichromaticMode(const mode::ParsBichromatic& p, AveragingConstructorParameters&&... a)
    : mode::Liouvillean<TEMPERATURE>(p.kappa,p.nTh),
      mode::Hamiltonian<true>(0,dcomp(mode::finiteTemperatureHamiltonianDecay(p,*this),-p.delta),p.eta,p.cutoff),
      ModeBase(p.cutoff,RF{"deltaOther",p.deltaOther,1},{CF{"(kappa*(2*nTh+1),delta)",conj(get_zI()),1},CF{"eta",get_eta(),sqrt(p.cutoff)},CF{"etaOther",p.etaOther,sqrt(p.cutoff)}},"Bichromatic mode"),
      AveragingType(std::forward<AveragingConstructorParameters>(a)...),
      zI_Other(mode::finiteTemperatureHamiltonianDecay(p,*this),-p.deltaOther)
  {
    mode::Hamiltonian<true>::getH_OverIs().push_back(
                                                     furnishWithFreqs(mode::pumping(p.etaOther,getDimension()),
                                                                      mode::mainDiagonal(zI_Other,getDimension()))
                                                     );
    getParsStream()<<"# Bichromatic pumping."; mode::isFiniteTempStream(getParsStream(),p.nTh,boost::mpl::bool_<TEMPERATURE>());
  }

private:
  const dcomp zI_Other;

};



template class BichromaticMode<true>;
template class BichromaticMode<false>;


#endif // CPPQEDELEMENTS_FREES_BICHROMATICMODE_H_INCLUDED
