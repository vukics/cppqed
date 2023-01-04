// Copyright Raimar Sandner 2012–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ParsParticleCavity_InterferenceWorkaround.h"

#include "Pars.h"

namespace particlecavity_interferenceworkaround {

ParsInterference::ParsInterference(parameters::Table& p, const std::string& mod)
  : modeInterference(p.addTitle("Interference",mod).add("modeInterference",mod,"Interference mode function",MFT_SIN)),
    kInterference(p.add<size_t>("kInterference",mod,"Interference mode function wave number",2)),
    uInterference(p.add("uInterference",mod,"Interference u parameter",1.)) {}

} // particlecavity_interferenceworkaround
