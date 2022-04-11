// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_FREES_HOMODYNEDMODE_H_INCLUDED
#define   CPPQEDELEMENTS_FREES_HOMODYNEDMODE_H_INCLUDED

#include "Mode_.h"
#include "TridiagonalHamiltonian.h"

#include "Pars.h"


namespace mode {

  
template<typename BASE>
struct ParsHomodyned : BASE // BASE either ParsLossy or ParsPumpedLossy
{
  typedef BASE BASE_class;

  dcomp &homodyneAmplitude;

  ParsHomodyned(parameters::Table& p, const std::string& mod="")
    : Pars(p,mod), BASE(p,mod),
      homodyneAmplitude(p.add("homodyneAmplitude",mod,"homodyneAmplitude",dcomp(0))) {}

};


class HomodynedBase
  : public Hamiltonian<false>, public structure::ElementLiouvillean<1,2>
{
protected:
  HomodynedBase(const ParsLossy& p, dcomp homodyneAmplitude, dcomp eta=0.);

  HomodynedBase(const ParsPumpedLossy& p, dcomp homodyneAmplitude) : HomodynedBase(p,homodyneAmplitude,p.eta) {}

private:
  void doActWithJ(NoTime, mode::StateVectorLow& psi, LindbladNo<0>) const;
  void doActWithJ(NoTime, mode::StateVectorLow& psi, LindbladNo<1>) const;

  double rate(NoTime, const mode::LazyDensityOperator&, LindbladNo<0>) const {return -1;}
  double rate(NoTime, const mode::LazyDensityOperator&, LindbladNo<1>) const {return -1;}

  const dcomp homodyneAmplitude_; ///< \f$\sqrt(f)\,e^{i\vartheta}\f$ in Charmichael’s notation

  const double kappa_, nTh_;

};


} // mode


template<typename AveragingType=mode::Averaged>
class HomodynedMode : public mode::HomodynedBase, public ModeBase, public AveragingType
{
public:
  template<typename BASE, typename... AveragingConstructorParameters>
  HomodynedMode(const mode::ParsHomodyned<BASE>& p, AveragingConstructorParameters&&... a)
    : mode::HomodynedBase(p,p.homodyneAmplitude),
      ModeBase(p.cutoff,{CF{"(kappa*(2*nTh+1),delta)",conj(get_zSch()),1},CF{"eta",get_eta(),sqrt(p.cutoff)},CF{"homodyneAmplitude",p.homodyneAmplitude,1.}},
               "Homodyned mode"),
      AveragingType(std::forward<AveragingConstructorParameters>(a)...)
  {
    getParsStream()<<"Homodyne detection.\n";
  }

};


namespace mode {


template<typename AveragingType, typename BASE, typename... AveragingConstructorParameters>
const Ptr make(const ParsHomodyned<BASE>& p, AveragingConstructorParameters&&... a)
{
  if (!isNonZero(p.homodyneAmplitude)) return make<AveragingType>(static_cast<const BASE&>(p),QMP_SCH,a...);
  else return std::make_shared<HomodynedMode<AveragingType> >(p,a...);
}

} // mode

#endif // CPPQEDELEMENTS_FREES_HOMODYNEDMODE_H_INCLUDED
