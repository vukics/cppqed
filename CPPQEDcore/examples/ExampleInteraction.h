// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ExampleMode.h"

#include "Interaction.h"
#include "TridiagonalHamiltonian.h"
#include "StateVector.h"

namespace basic {

class InteractionX_X
  : public structure::Interaction<2>, public TridiagonalHamiltonian<2,false>
{
public:
  InteractionX_X(std::shared_ptr<PumpedLossyMode>, std::shared_ptr<PumpedLossyMode>, double g);

};

} // basic


namespace hierarchical {

class InteractionX_X
  : public structure::Interaction<2>, public TridiagonalHamiltonian<2,true>
{
public:
  InteractionX_X(ModeBase::Ptr, ModeBase::Ptr, double g);

};

} // hierarchical


namespace hierarchical {

class InteractionX_X_Correlations : public InteractionX_X, public structure::ElementAveraged<2>
{
public:
  InteractionX_X_Correlations(ModeBase::Ptr m0, ModeBase::Ptr m1, double g)
    : InteractionX_X(m0,m1,g),
      ElementAveraged<2>("ModeCorrelations",{"<XQ>","<XP>","<YQ>","<YP>"})
  {}

private:
  const Averages average_v(NoTime, const quantumdata::LazyDensityOperator<2>&) const;

};

} // hierarchical
