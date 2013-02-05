// -*- C++ -*-
#ifndef STRUCTURE_DYNAMICSBASE_H_INCLUDED
#define STRUCTURE_DYNAMICSBASE_H_INCLUDED

#include "DynamicsBaseFwd.h"

#include "ComplexExtensions.h"

#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>
#include <boost/tuple/tuple.hpp>

#include <boost/assign/list_of.hpp>

#include <list>


#define TTD_FREQUENCY_MAP(T) std::list<boost::tuple<std::string,T,double> >


#define FREQS boost::assign::tuple_list_of


namespace structure {


class DynamicsBase : private boost::noncopyable
{
public:
  typedef boost::shared_ptr<const DynamicsBase> Ptr;

  typedef TTD_FREQUENCY_MAP(double)    RealFreqs;
  typedef TTD_FREQUENCY_MAP(dcomp ) ComplexFreqs;

  explicit DynamicsBase(const    RealFreqs& =RealFreqs   (), 
                        const ComplexFreqs& =ComplexFreqs());

  // Calculating the fastest timescale of the system
  double highestFrequency() const;

  // Displaying parameters of the system
  std::ostream& displayParameters(std::ostream&) const;

  virtual ~DynamicsBase() {}

protected:
  std::stringstream& getParsStream() {return paramsStream_;}
  
  RealFreqs&    getRealFreqs   () {return    realFreqs_;}
  ComplexFreqs& getComplexFreqs() {return complexFreqs_;}

  virtual std::ostream& displayMoreParameters(std::ostream&) const;

private:
  RealFreqs       realFreqs_;
  ComplexFreqs complexFreqs_;

  std::stringstream paramsStream_;

};


} // structure


#undef TTD_FREQUENCY_MAP


#endif // STRUCTURE_DYNAMICSBASE_H_INCLUDED
