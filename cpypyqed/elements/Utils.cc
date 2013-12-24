// -*- C++ -*-

#include "PythonExtension.h"
#include "Namespaces.h"

#include "ModeFunction.h"
#include "ParticleInitialCondition.h"
#include "QM_PictureFwd.h"

using namespace boost::python;

namespace pythonext{

void export_utils()
{
  enum_<ModeFunctionType>("ModeFunctionType")
    .value("COS",  MFT_COS)
    .value("SIN",  MFT_SIN)
    .value("PLUS", MFT_PLUS)
    .value("MINUS",MFT_MINUS)
  ;
  enum_<QM_Picture>("QMP")
    .value("IP", QMP_IP)
    .value("UIP", QMP_UIP)
    .value("SCH", QMP_SCH)
  ;
  {
    scope namespaceScope = particleNameSpace;
    class_<particle::InitialCondition>("InitialCondition", init<double, double, double, optional<bool> >())
      .add_property("getX0", &particle::InitialCondition::getX0)
      .add_property("getK0", &particle::InitialCondition::getK0)
      .add_property("getSig", &particle::InitialCondition::getSig)
      .add_property("isInK", &particle::InitialCondition::isInK)
    ;
  }
}

} // pythonext