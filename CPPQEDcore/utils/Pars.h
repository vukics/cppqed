// -*- C++ -*-
#ifndef UTILS_PARS_H_INCLUDED
#define UTILS_PARS_H_INCLUDED

#include "ParsFwd.h"

#include "BooleanNegatedProxy.h"
#include "Exception.h"

#include <boost/ptr_container/ptr_list.hpp>

#include <string>
#include <sstream>


namespace parameters {

const size_t maxTypeLabelLength=24;


void update(ParameterTable&, int, char*[], const std::string& ="--");


//////////////////
// 
// Pars Exceptions
//
//////////////////


class Exception : public cpputils::Exception
{
public:
  virtual ~Exception() {}
};


class NamedException : public Exception 
{
public:
  NamedException(const std::string& name);

  const std::string& getName() const {return name_;}

private:
  const std::string name_;
};


class UnrecognisedParameterException : public NamedException 
{
public: 
  UnrecognisedParameterException(const std::string& name) : NamedException(name) {}
};


class AttemptedRecreationOfParameterException : public NamedException 
{
public: 
  AttemptedRecreationOfParameterException(const std::string& name);
};



////////////////////
//
// Auxiliary Classes
//
////////////////////

// Common Base Class

class ParameterBase
{
public:
  ParameterBase(const std::string& s, const std::string& d) : s_(s), d_(d) {}
  virtual ~ParameterBase() {}

  const std::string& getS() const {return s_;}
  const std::string& getD() const {return d_;}
    
  virtual void print(size_t, size_t, size_t) const = 0;
  virtual void read (std::istream&)                = 0;

private:
  const std::string s_; // as appearing in command line
  const std::string d_; // short description

};


// Template containing value for the Parameter

template<typename T>
class Parameter : public ParameterBase 
{
public:
  Parameter(const std::string& s, const std::string& d, const T& v) : ParameterBase(s,d), v_(v) {}
  
  void print(size_t smw, size_t tmw, size_t dmw) const;
  void read(std::istream& is) {is>>v_;}

  const T& getReference() const {return v_;}

  T& getReference() {return const_cast<T&>(static_cast<const Parameter*>(this)->getReference());}

private:
  T v_;
  
};


//////////////////
//
// Parameter Table
//
//////////////////

class ParameterTable
{
public:
  typedef boost::ptr_list<ParameterBase> Impl;

  ParameterTable() : table_(), smwidth_(0), tmwidth_(6), dmwidth_(0), stream_() {} // tmwidth_ cf bool!
  ~ParameterTable() {}

  const ParameterBase& operator[](const std::string&) const;

  ParameterBase& operator[](const std::string& s) {return const_cast<ParameterBase&>(static_cast<const ParameterTable&>(*this)[s]);}

  template<typename T> T& add(const std::string& s, const std::string& d, const T& v);

  bool& add(const std::string& s, const std::string& d, bool v);

  template<typename T> T& addMod(const std::string& s, const std::string& mod, const std::string& d, const T& v)
  {
    return add(s+mod,d,v);
  }

  ParameterTable& addTitle(const std::string& s, const std::string& mod="");

  void printList() const;
  
  std::iostream& getStream() {return stream_;}

private:
  Impl table_;

  size_t smwidth_; // maximal width of s_ entries
  size_t tmwidth_; // maximal width of typeIDs
  size_t dmwidth_; // maximal width of d_ entries

  std::stringstream stream_;
  
};


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



} // parameters


#endif // UTILS_PARS_H_INCLUDED
