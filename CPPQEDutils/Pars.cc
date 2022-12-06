// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Pars.h"

// #include "BooleanNegatedProxy.h"

#include "Version.h"


using namespace popl;


std::shared_ptr<Switch> help_switch, version_switch;


OptionParser optionParser(std::string pre, std::string post)
{
  OptionParser res{pre+versionHelper()+post};
  help_switch=res.add<Switch>("h", "help", "list physical parameters & program options");
  version_switch=res.add<Switch>("v", "version", "display version information");
  return res;
}

//    if (temp=="help") {
//      table.printList(); exit(0);
//    }
//    if (temp=="version") {
//      cerr << versionHelper(); exit(0);
//    }
// }

/*

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
}*/
