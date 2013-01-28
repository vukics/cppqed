// -*- C++ -*-
#ifndef   UTILS_INCLUDE_PARSFWD_H_INCLUDED
#define   UTILS_INCLUDE_PARSFWD_H_INCLUDED

#include<string>


namespace parameters {

class Exception;
class NamedException;
class UnrecognisedParameterException;
class AttemptedRecreationOfParameterException;

template<typename>
class Parameter;

class ParameterTable;

struct TitleLine {}; 
// A tagging class for introducing dummy parameters into the Table, which simply create a newline and a title at the listing.


} // parameters


#endif // UTILS_INCLUDE_PARSFWD_H_INCLUDED
