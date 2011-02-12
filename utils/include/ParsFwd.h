// -*- C++ -*-
#ifndef   _PARS_FWD_H
#define   _PARS_FWD_H

#include<string>


namespace parameters {

class ParsException;
class ParsNamedException;
class UnrecognisedParameterException;
class AttemptedRecreationOfParameterException;
class ParameterTypeMismatchException;

template<typename>
class Parameter;

class ParameterTable;

struct TitleLine {}; 
// A tagging class for introducing dummy parameters into the Table, which simply create a newline and a title at the listing.


} // parameters


#endif // _PARS_FWD_H
