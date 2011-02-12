// -*- C++ -*-
#ifndef _MODE_FUNCTION_TYPE_H
#define _MODE_FUNCTION_TYPE_H

#include "ModeFunctionTypeFwd.h"

#include "ComplexExtensions.h"

#include<iosfwd>

inline bool isComplex(ModeFunctionType mf) {return (mf==MFT_PLUS || mf==MFT_MINUS);}

std::ostream& operator<<(std::ostream&, ModeFunctionType);
std::istream& operator>>(std::istream&, ModeFunctionType&);


const dcomp modeFunction(ModeFunctionType,double);


#endif // _MODE_FUNCTION_TYPE_H
