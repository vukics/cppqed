// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "KeyPrinter.h"

#include <iostream>
#include <iomanip>

using namespace std;


ostream& cppqedutils::KeyPrinter::stream(ostream& os, size_t& i) const
{
  os<<keyTitle_;
  for (const auto& kl : keyLabels_) os<<endl<<setw(2)<<i++<<". "<<kl;
  return os<<endl;
}

