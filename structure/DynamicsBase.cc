#include "DynamicsBase.h"

#include "impl/FormDouble.tcc"
#include "Range.h"

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <boost/mpl/identity.hpp>

namespace mpl=boost::mpl;

using namespace std;

namespace structure {


DynamicsBase::DynamicsBase(const RealFreqs& realFreqs, const ComplexFreqs& complexFreqs) 
  : realFreqs_(realFreqs), complexFreqs_(complexFreqs), paramsStream_(std::stringstream::out)
{
  paramsStream_.precision(formdouble::actualPrecision(FormDouble::overallPrecision));
}



#define TTD_NAMED_FREQUENCY(T) std::list<boost::tuple<std::string,T,double> >::value_type


namespace {

inline double absReal   (const TTD_NAMED_FREQUENCY(double)& p) {return fabs(p.get<1>()*p.get<2>());}
inline double absComplex(const TTD_NAMED_FREQUENCY(dcomp )& p) {return  abs(p.get<1>()*p.get<2>());}


template<typename T>
bool templateCompare(const typename TTD_NAMED_FREQUENCY(T)& p1, const typename TTD_NAMED_FREQUENCY(T)& p2,
		     boost::function<double(const typename TTD_NAMED_FREQUENCY(T)&)> f)
{
  return f(p1)<f(p2);
}

} // unnamed namespace


double DynamicsBase::highestFrequency() const
{
  using boost::max_element;
  return
    max(
	realFreqs_.size() 
	? 
	absReal   (*max_element(realFreqs_   ,bind(templateCompare<double>,_1,_2,absReal   )))
	:
	0,
	complexFreqs_.size()
	?
	absComplex(*max_element(complexFreqs_,bind(templateCompare<dcomp >,_1,_2,absComplex)))
	:
	0
	);

}



void DynamicsBase::displayParameters(ostream& os) const
{
  os<<paramsStream_.str();
  displayMoreParameters(os);
  os<<endl;
}


namespace {

template<typename T>
void displayFreq(ostream& os, int precision, const typename TTD_NAMED_FREQUENCY(T)& pair)
{
  os<<"# "<<pair.template get<0>()<<"="<<formdouble::zeroWidth(precision)(pair.template get<1>())<<endl;
}

} // unnamed namespace


void DynamicsBase::displayMoreParameters(ostream& os) const
{
  using boost::ref; using boost::for_each;
  for_each(   realFreqs_,bind(displayFreq<double>,ref(os),os.precision(),_1));
  for_each(complexFreqs_,bind(displayFreq<dcomp >,ref(os),os.precision(),_1));  
}


} // structure


#undef TTD_NAMED_FREQUENCY
