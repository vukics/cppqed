// Copyright András Vukics 2006–2015. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef CPPQEDELEMENTS_UTILS_MODEFUNCTION_H_INCLUDED
#define CPPQEDELEMENTS_UTILS_MODEFUNCTION_H_INCLUDED

#include "ComplexExtensions.h"

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

#include <iosfwd>

enum ModeFunctionType {MFT_SIN, MFT_COS, MFT_PLUS, MFT_MINUS};

inline bool isComplex(ModeFunctionType mf) {return (mf==MFT_PLUS || mf==MFT_MINUS);}

std::ostream& operator<<(std::ostream&, ModeFunctionType);
std::istream& operator>>(std::istream&, ModeFunctionType&);


const dcomp modeFunction(ModeFunctionType,double);


typedef boost::tuple<ModeFunctionType,ptrdiff_t> ModeFunction;

// std::ostream& operator<<(std::ostream&, const particle::ModeFunction&);

#endif // CPPQEDELEMENTS_UTILS_MODEFUNCTION_H_INCLUDED
