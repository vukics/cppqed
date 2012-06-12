// -*- C++ -*-
#ifndef _STRUCTURE_DYNAMICS_BASE_H
#define _STRUCTURE_DYNAMICS_BASE_H

#include "DynamicsBaseFwd.h"

#include "ComplexExtensions.h"

#include<boost/utility.hpp>

#include<boost/tuple/tuple.hpp>

#include<list>


#define TTD_FREQUENCY_MAP(T) std::list<boost::tuple<std::string,T,double> >



namespace structure {


class DynamicsBase : private boost::noncopyable
{
public:
  typedef TTD_FREQUENCY_MAP(double)    RealFreqs;
  typedef TTD_FREQUENCY_MAP(dcomp ) ComplexFreqs;

  explicit DynamicsBase(const    RealFreqs& =RealFreqs   (), 
			const ComplexFreqs& =ComplexFreqs());

  // Calculating the fastest timescale of the system
  double highestFrequency() const;

  // Displaying parameters of the system
  void displayParameters(std::ostream&) const;

  virtual ~DynamicsBase() {}

protected:
  std::stringstream& getParsStream() {return paramsStream_;}
  
  RealFreqs&    getRealFreqs   () {return    realFreqs_;}
  ComplexFreqs& getComplexFreqs() {return complexFreqs_;}

  virtual void displayMoreParameters(std::ostream&) const;

private:
  RealFreqs       realFreqs_;
  ComplexFreqs complexFreqs_;

  std::stringstream paramsStream_;

};


} // structure


#undef TTD_FREQUENCY_MAP


#endif // _STRUCTURE_DYNAMICS_BASE_H
