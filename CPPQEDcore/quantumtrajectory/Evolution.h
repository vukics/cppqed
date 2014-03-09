/// \briefFile{Convenience header collecting to one place most of the components of C++QED core needed for writing scripts (basically, the core part of the Level-1 interface)}
#ifndef CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION_H_INCLUDED
#define CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION_H_INCLUDED

#include "BlitzArrayTraits.h"

#include "Evolution.tcc"
#include "QM_Picture.h"
#include "Tridiagonal.tcc"

#include "EvolvedGSL.tcc"
#include "Pars.tcc"

#include "component_versions.h"


class ParameterTable : public parameters::ParameterTable {};


void update(ParameterTable& p, int argc, char* argv[], const std::string& prefix="--")
{
  updateVersionstring(cppqed_component_versions());
  parameters::update(p,argc,argv,prefix);
}


/// Convenience version of parameters::update meant to tackle the problem described in Sec. \ref masterequationlimitations
QM_Picture& updateWithPicture(ParameterTable& p, int argc, char* argv[], const std::string& prefix="--")
{
  QM_Picture& qmp=p.add("picture","Quantum mechanical picture",QMP_IP);
  update(p,argc,argv,prefix);
  try {
    const evolution::Method method=dynamic_cast<const parameters::Parameter<evolution::Method>&>(p["evol"]).get();
    if ((method==evolution::MASTER || method==evolution::MASTER_FAST) && qmp==QMP_IP) qmp=QMP_UIP;
  } catch (const parameters::UnrecognisedParameterException&) {}
  return qmp;
}


#endif // CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION_H_INCLUDED

