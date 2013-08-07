#include "KeyPrinter.h"

#include "Range.h"

#include <boost/lambda/lambda.hpp>

#include <iostream>
#include <iomanip>

using namespace std;



cpputils::KeyPrinter::KeyPrinter(const string& keyTitle, const KeyLabels& keyLabels) 
  : keyTitle_(keyTitle), keyLabels_(keyLabels) 
{
}



ostream& cpputils::KeyPrinter::displayKey(ostream& os, size_t& i) const
{
  namespace bll=boost::lambda;

  os<<"# "<<keyTitle_;
  boost::for_each(keyLabels_,os<<bll::constant("\n# ")<<bll::constant(setw(2))<<bll::var(i)++<<". "<<bll::_1);
  return os<<endl;
}

