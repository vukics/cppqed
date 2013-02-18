#include "impl/Pars.tcc"

#include "Range.h"

#include "BooleanNegatedProxy.h"

#include <blitz/bzconfig.h>

#include <gsl/gsl_version.h>

#include <boost/bind.hpp>
#include <boost/version.hpp>

#include <boost/preprocessor/stringize.hpp>

#include <iostream>
#include <iomanip>
#include <fstream>



using namespace std;
using boost::cref;


parameters::NamedException::NamedException(const std::string& name)
  : name_(name)
{}


parameters::AttemptedRecreationOfParameterException::AttemptedRecreationOfParameterException(const std::string& name)
  : parameters::NamedException(name)
{
  cerr<<name<<endl;
}


namespace parameters {


namespace {

const string errorMessage(const string& sc)
{
  return "\nSee "+sc+"help for a list of parameters\nSupply all parameters with "+sc+" prefix\n\n";
}

void versionHelper()
{
  string versionString="C++QED Version 2 Milestone 9";
  {
    ifstream file((string(BOOST_PP_STRINGIZE(TOP_LEVEL_DIRECTORY))+string("/.bzr/branch/last-revision")).c_str());
    if (file.is_open()) {
      versionString="Revision #";
      char temp[256];
      file.getline(temp,256);
      versionString+=string(temp);
    }
  }

  cout<<endl<<versionString<<"\nAndras Vukics, vukics@users.sourceforge.net\n\nCompiled with\nBoost library collection : Version "<<string(BOOST_LIB_VERSION).replace(1,1,".")<<"\nGnu Scientific Library   : Version "<<GSL_VERSION<<"\nBlitz++ numerical library: Config date: "<<BZ__config_date<<endl<<endl;
}

}


const ParameterBase&
ParameterTable::operator[](const string& s) const
{
  Impl::const_iterator i=boost::find_if(table_,bind(equal_to<string>(),bind(&ParameterBase::getS,_1),cref(s)));

  if (i==table_.end()) throw UnrecognisedParameterException(s);
  return *i;
}


void ParameterTable::printList() const
{
  cout<<"\nhelp       display this list\nversion    display version information\n\n";
  boost::for_each(table_,bind(&ParameterBase::print,_1,smwidth_,tmwidth_,dmwidth_));
  cout<<endl;
}


bool& ParameterTable::add(const std::string& s, const std::string& d, bool v)
{
  bool& res=add<bool>(s,d,v);
  add("no_"+s,d,cpputils::BooleanNegatedProxy(res));
  return res;
}


ParameterTable& ParameterTable::addTitle(const std::string& s, const std::string& mod)
{
  try {(*this)[s+mod]; throw AttemptedRecreationOfParameterException(s);}
  catch (UnrecognisedParameterException) {
    table_.push_back(new Parameter<TitleLine>(s+mod,"",TitleLine()));
  }
  return *this;
}


void update(ParameterTable& table, int argc, char* argv[], const string& sc)
{
  if (argc<2) return;

  iostream& line=table.getStream();
  for (int i=1; i<argc; ++i) line<<' '<<argv[i];

  string temp;
  size_t scl=sc.length();
  while (line>>temp) {
    if (string(temp,0,scl)!=sc) {
      cerr<<"\nSyntax error in command line around \""<<temp<<"\""<<errorMessage(sc);
      abort();
    }
    temp=string(temp,scl,temp.length()-scl);
    if (temp=="help") {
      table.printList(); exit(0);
    }
    if (temp=="version") {
      versionHelper(); exit(0);
    }
    else
      try {
	table[temp].read(line);
	if (line.fail()) {
	  cerr<<"\nParameter \""<<temp<<"\" supplied with incorrect syntax"<<errorMessage(sc);
	  abort();
	}
      }
      catch (const UnrecognisedParameterException& urp) {
	cerr<<"\nProblem in command line around \""<<urp.getName()<<'\"'<<errorMessage(sc);
	abort();
      }
  }
}






template<>
void Parameter<bool>::print(size_t smw, size_t tmw, size_t dmw) const
{
  using namespace std;
  cout<<setw(smw+3)<<left<<getS()
      <<setw(tmw+3)<<left<<"switch"
      <<setw(dmw+3)<<left<<getD()<<v_<<endl;
}


template<>
void Parameter<bool>::read(std::istream&)
{
  v_=true;
}


template<>
void Parameter<cpputils::BooleanNegatedProxy>::read(std::istream&)
{
  v_=true;
}


template<>
void Parameter<TitleLine>::print(size_t, size_t, size_t) const
{
  using namespace std;
  cout<<endl<<"*** "<<getS()<<endl;
}


template<>
void Parameter<TitleLine>::read(std::istream&)
{
  throw UnrecognisedParameterException(getS());
}


} // parameters
