// -*- C++ -*-

#include "PythonExtension.h"
#include "Namespaces.h"

#include "ModeFunction.h"
#include "ParticleInitialCondition.h"
#include "QM_Picture.h"

#include <boost/python/docstring_options.hpp>

using namespace boost::python;

namespace pythonext{

QM_Picture updateWithPicture(parameters::ParameterTable &p, list args, str s=str("--")){
  int argc = len(args);
  char **argv = new char*[argc];
  for (int c=0; c<argc; c++) {
    argv[c] = extract<char *>(args[c]);
  }
  QM_Picture result=picture::updateWithPicture(p,argc,argv,extract<std::string>(s));
  delete[] argv;
  return result;
}

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
  def("isComplex", &isComplex);
  {
    scope namespaceScope = particleNameSpace;
    class_<particle::InitialCondition>("InitialCondition", init<double, double, double, optional<bool> >())
      .def("getX0", &particle::InitialCondition::getX0)
      .def("getK0", &particle::InitialCondition::getK0)
      .def("getSig", &particle::InitialCondition::getSig)
      .def("isInK", &particle::InitialCondition::isInK)
    ;
  }
  {
    docstring_options docOpts;
    docOpts.disable_cpp_signatures();
    def("updateWithPicture", &updateWithPicture);
    def("updateWithPicture", &updateWithPicture,
R"doc(Wrapper of :elements:`picture::updateWithPicture`. Note that the signature differs a little bit
compared to the C++ version, it is not necessary to pass in the argument list length.

:param p: The :class:`ParameterTable`
:type p: :class:`ParameterTable`
:param list argv: The argument list
:param str prefix: The prefix (default "``--``")
:returns qmp: The quantum mechanical picture to use)doc",
        (arg("p"),"argv",arg("prefix")="--")
    );
  }
}

} // pythonext