// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Pars.h"

// #include "BooleanNegatedProxy.h"

#include "Version.h"


using namespace popl;


static std::shared_ptr<Switch> help_switch, version_switch;


OptionParser optionParser(std::string pre, std::string post)
{
  OptionParser res{pre+versionHelper()+post+"Command-line arguments"};
  help_switch=res.add<Switch>("h", "help", "list physical parameters & program options");
  version_switch=res.add<Switch>("v", "version", "display version information");
  return res;
}


void parse(popl::OptionParser& op, int argc, const char* const argv[])
{
  for (const char* const* i=argv;  i<argv+argc; ++i) (parsedCommandLine+=*i)+=' ';

  if ( !bool(help_switch) || !bool(version_switch) )
    throw std::logic_error("Uninitialized help and/or version switches – did you maybe forget using the optionParser() function?");
 
  op.parse(argc, argv);
  
  using std::cerr;

  // print auto-generated help message
  if (version_switch->count()) {
    cerr << versionHelper(); exit(0);
  }
  
  // print auto-generated help message
  if (help_switch->count()) {
    if (help_switch->count() == 1)
      cerr << op << "\n";
    else if (help_switch->count() == 2)
      cerr << op.help(Attribute::advanced) << "\n";
    else if (help_switch->count() > 2)
      cerr << op.help(Attribute::expert) << "\n";
    exit(0);
  }

  // show all non option arguments (those without "-o" or "--option")
  for (const auto& non_option_arg: op.non_option_args())
    cerr << "non_option_args: " << non_option_arg << "\n";

  // show unknown options (undefined ones, like "-u" or "--undefined")
  for (const auto& unknown_option: op.unknown_options())
    cerr << "unknown_options: " << unknown_option << "\n";
}
