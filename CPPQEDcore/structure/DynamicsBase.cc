// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "DynamicsBase.h"
#include "FormDouble.h"


using namespace std;

namespace structure {


DynamicsBase::DynamicsBase(const RealFreqs& realFreqs, const ComplexFreqs& complexFreqs) 
  : realFreqs_(realFreqs), complexFreqs_(complexFreqs), paramsStream_()
{
  paramsStream_.precision(formdouble::actualPrecision(FormDouble::overallPrecision));
}



template <typename T>
inline double absFreq(const std::tuple<std::string,T,double>& p) {return std::abs(get<1>(p)*get<2>(p));}


double DynamicsBase::highestFrequency() const
{
  const auto compareFreqs=[](const auto& a, const auto& b) {return absFreq(a)<absFreq(b);};
  
  return max(
    realFreqs_.size() ? absFreq(*std::ranges::max_element(realFreqs_,compareFreqs)) : 0,
    complexFreqs_.size() ? absFreq(*std::ranges::max_element(complexFreqs_,compareFreqs)) : 0);

}



std::ostream& DynamicsBase::streamParameters(ostream& os) const
{
  return streamMoreParameters(os<<paramsStream_.str())<<endl;
}


std::ostream& DynamicsBase::streamMoreParameters(ostream& os) const
{
  for(const auto& rf :    realFreqs_) os<<get<0>(rf)<<"="<<formdouble::zeroAdditional(os.precision())(get<1>(rf))<<endl;
  for(const auto& cf : complexFreqs_) os<<get<0>(cf)<<"="<<formdouble::zeroAdditional(os.precision())(get<1>(cf))<<endl;
  return os;
}


} // structure

