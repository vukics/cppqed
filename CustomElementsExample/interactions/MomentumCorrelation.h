#ifndef INTERACTIONS_MOMENTUMCORRELATION_H_INCLUDED
#define INTERACTIONS_MOMENTUMCORRELATION_H_INCLUDED

#include "ElementAveraged.h"
#include "Interaction.h"
#include "Particle.h"
#include "LazyDensityOperator.h"

class MomentumCorrelation : public structure::ElementAveraged<2>
{
protected:
  typedef structure::ElementAveraged<2> EA_Base;
  typedef quantumdata::LazyDensityOperator<2> LazyDensityOperator;
  MomentumCorrelation();
  
private:
  const Averages average_v(const LazyDensityOperator&) const;
  void           process_v(Averages&)                  const {}
};

namespace momentumcorrelationinteraction {

class Base : public structure::Interaction<2>
{
protected:
  typedef structure::Interaction<2> IA_Base;
  
  Base(particle::Ptr p1, particle::Ptr p2);
};

} // momentumcorrelationinteraction

#define BIG_NAMESPACE_NAME             momentumcorrelationinteraction
#define BIG_CLASS_NAME                 MomentumCorrelationInteraction
#define BIG_ADDITIONAL_PARAMETERS      
#define BIG_ADDITIONAL_PARAMETERS_PASS 

#include "details_BinaryInteractionGenerator.h"


#endif // INTERACTIONS_MOMENTUMCORRELATION_H_INCLUDED
