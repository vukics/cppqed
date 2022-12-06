// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "popl.hpp"


popl::OptionParser optionParser(std::string pre="", std::string post="");

void parse(const popl::OptionParser&);
