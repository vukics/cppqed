// -*- C++ -*-
#ifndef   UTILS_PARS_TCC_INCLUDED
#define   UTILS_PARS_TCC_INCLUDED

#include "Pars.h"

#include <iostream>
#include <iomanip>
#include <typeinfo>


namespace parameters {


template<typename T>
void
Parameter<T>::print_v(size_t smw, size_t tmw, size_t dmw) const
{
  using namespace std;
  const string name(typeid(T).name());
  cout<<setw(smw+3)<<left<<getS()
      <<setw(tmw+3)<<left<<(name.size()>maxTypeLabelLength-3 ? name.substr(0,maxTypeLabelLength-3)+"..." : name)
      <<setw(dmw+3)<<left<<getD()<<v_<<endl;
}


template<typename T>
T&
ParameterTable::add(const std::string& s, const std::string& d, const T& v)
{
  using namespace std;
  try {(*this)[s]; throw AttemptedRecreationOfParameterException(s);}
  catch (UnrecognisedParameterException) {
    Parameter<T>* pptr=new Parameter<T>(s,d,v);
    table_.push_back(pptr);
    smwidth_=max(smwidth_,s.length());
    tmwidth_=max(tmwidth_,min(strlen(typeid(T).name()),maxTypeLabelLength));
    dmwidth_=max(dmwidth_,d.length());
    return pptr->getReference();
  }
}


} // parameters

#endif // UTILS_PARS_TCC_INCLUDED
