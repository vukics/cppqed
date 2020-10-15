// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Pars.h"

#include "BooleanNegatedProxy.h"

#include "IO_Manip.h"
#include "Version.h"

#include <boost/range/algorithm.hpp>

#include <boost/bind.hpp>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>


using namespace std;


const parameters::Base&
parameters::Table::operator[](const string& s) const
{
  auto i=boost::find_if(table_,bind(equal_to<string>(),bind(&Base::getS,_1),boost::cref(s)));

  if (i==table_.end()) throw UnrecognisedParameterException(s);
  return *i;
}


parameters::Table::Table() : table_(), smwidth_(0), tmwidth_(6), dmwidth_(0), parsedCommandLine_(std::make_shared<std::string>("")) // tmwidth_ cf bool!
{
  IO_Manipulator::_(cout);
  IO_Manipulator::_(cerr);
} 


void parameters::Table::printList() const
{
  cout<<"\nhelp       display this list\nversion    display version information\n\n";
  boost::for_each(table_,bind(&Base::print,_1,smwidth_,tmwidth_,dmwidth_));
  cout<<endl;
}


///////////////////////////////////
//
// Boolean Typed Specialization
//
///////////////////////////////////


template<>
void parameters::Typed<bool>::print_v(size_t smw, size_t tmw, size_t dmw) const
{
  using namespace std;
  cout<<setw(smw+3)<<left<<getS()
      <<setw(tmw+3)<<left<<"switch"
      <<setw(dmw+3)<<left<<getD()<<v_<<endl;
}


template<>
void parameters::Typed<bool>::read_v(std::istream&)
{
  v_=true;
}


template<>
void parameters::Typed<cpputils::BooleanNegatedProxy>::read_v(std::istream&)
{
  v_=true;
}


bool& parameters::Table::add(const std::string& s, const std::string& d, bool v)
{
  bool& res=add<bool>(s,d,v);
  add("no_"+s,d,cpputils::BooleanNegatedProxy(res));
  return res;
}


////////////////////////////
//
// TitleLine specializations
//
////////////////////////////

namespace {

// A tagging class for introducing dummy parameters into the Table, which simply create a newline and a title at the listing.
struct TitleLine {}; 

}

template<>
void parameters::Typed<TitleLine>::print_v(size_t, size_t, size_t) const
{
  using namespace std;
  cout<<endl<<"*** "<<getS()<<endl;
}


template<>
void parameters::Typed<TitleLine>::read_v(std::istream&)
{
  throw UnrecognisedParameterException(getS());
}


parameters::Table& parameters::Table::addTitle(const std::string& s, const std::string& mod)
{
  try {(*this)[s+mod]; throw AttemptedRecreationOfParameterException(s);}
  catch (const UnrecognisedParameterException&) {
    table_.push_back(new Typed<TitleLine>(s+mod,"",TitleLine()));
  }
  return *this;
}


void parameters::update(parameters::Table& table, int argc, char* argv[], const string& prefix)
{
  struct ErrorMessage {
    static const string _(const string& prefix) {return "\nSee "+prefix+"help for a list of parameters\nSupply all parameters with "+prefix+" prefix\n\n";}
  };

  std::stringstream line; IO_Manipulator::_(line);
  for (char** i=argv+1;  i<argv+argc; ++i) line<<' '<<*i;
  table.setParsedCommandLine(string(*argv)+line.str());

  if (argc<2) return;

  string temp;
  size_t prefixl=prefix.length();
  while (line>>temp) {
    if (string(temp,0,prefixl)!=prefix) {
      cerr<<"\nSyntax error in command line around \""<<temp<<"\""<<ErrorMessage::_(prefix);
      abort();
    }
    temp=string(temp,prefixl,temp.length()-prefixl);
    if (temp=="help") {
      table.printList(); exit(0);
    }
    if (temp=="version") {
      cerr << versionHelper(); exit(0);
    }
    else
      try {
        table[temp].read(line);
        if (line.fail()) {
          cerr<<"\nParameter \""<<temp<<"\" supplied with incorrect syntax"<<ErrorMessage::_(prefix);
          abort();
        }
      }
      catch (const UnrecognisedParameterException& urp) {
        cerr<<"\nProblem in command line around \""<<urp.what()<<'\"'<<ErrorMessage::_(prefix);
        abort();
      }
  }
}
