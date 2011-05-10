// -*- C++ -*-
#ifndef _MODE_FUNCTION_TYPE_H
#define _MODE_FUNCTION_TYPE_H

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

#endif // _MODE_FUNCTION_TYPE_H
