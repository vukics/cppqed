#include "ParsParticleCavity_InterferenceWorkaround.h"

#include "Pars.h"

namespace particlecavity_interferenceworkaround {

ParsInterference::ParsInterference(parameters::ParameterTable& p, const std::string& mod)
  : modeInterference(p.addTitle("Interference",mod).addMod("modeInterference",mod,"Interference mode function",MFT_SIN)),
    kInterference(p.addMod<size_t>("kInterference",mod,"Interference mode function wave number",2)),
    uInterference(p.addMod("uInterference",mod,"Interference u parameter",1.)) {}

} // particlecavity_interferenceworkaround
