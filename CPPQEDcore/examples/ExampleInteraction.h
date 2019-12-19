// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ExampleMode.h"

#include "Interaction.h"
#include "TridiagonalHamiltonian.tcc"
#include "StateVector.tcc"

namespace basic {

class InteractionX_X
  : public Interaction<2>, public TridiagonalHamiltonian<2,false>
{
public:
  InteractionX_X(const PumpedLossyMode&, const PumpedLossyMode&, double g);

};

} // basic


namespace hierarchical {

class InteractionX_X
  : public Interaction<2>, public TridiagonalHamiltonian<2,true>
{
public:
  InteractionX_X(const ModeBase&, const ModeBase&, double g);

};

} // hierarchical


namespace hierarchical {

class InteractionX_X_Correlations : public InteractionX_X, public ElementAveraged<2>
{
public:
  InteractionX_X_Correlations(const ModeBase& m0, const ModeBase& m1, double g)
    : InteractionX_X(m0,m1,g),
      ElementAveraged<2>("ModeCorrelations",{"<XQ>","<XP>","<YQ>","<YP>"})
  {}

private:
  const Averages average_v(NoTime, const LazyDensityOperator&) const;

};

} // hierarchical