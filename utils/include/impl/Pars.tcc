// -*- C++ -*-
#ifndef   _PARS_IMPL_H
#define   _PARS_IMPL_H

#include<iostream>
#include<iomanip>

#include "BooleanNegatedProxy.h"

namespace parameters {


template<typename T>
void
Parameter<T>::print(size_t smw, size_t tmw, size_t dmw) const
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


///////////////////////////////////
//
// Boolean Parameter Specialization
//
///////////////////////////////////


template<>
void Parameter<bool>::print(size_t smw, size_t tmw, size_t dmw) const;


template<>
void Parameter<bool>::read(std::istream&);


template<>
void Parameter<cpputils::BooleanNegatedProxy>::read(std::istream&);


////////////////////////////
//
// TitleLine specializations
//
////////////////////////////


template<>
void Parameter<TitleLine>::print(size_t smw, size_t tmw, size_t dmw) const;


template<>
void Parameter<TitleLine>::read(std::istream&);


}

#endif // _PARS_IMPL_H
