// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "popl.hpp"


popl::OptionParser optionParser(std::string pre="", std::string post="");

void parse(const popl::OptionParser&);

template <typename T>
auto& add(popl::OptionParser& op, std::string option, std::string description, T defaultValue, auto&& binding)
{
  op.add<popl::Value<T>>("",option,description,defaultValue,std::forward<decltype(binding)>(binding));
  return op;
}

template <typename T>
auto& add(popl::OptionParser& op, std::string option, std::string mod, std::string description, T defaultValue, auto&& binding)
{
  return add(op,option+mod,description,defaultValue,std::forward<decltype(binding)>(binding));
}

/// A no-op. TODO: resurrect the add-title feature
auto& addTitle(popl::OptionParser& op, std::string /* title */, std::string /* mod */ = "") {return op;}
