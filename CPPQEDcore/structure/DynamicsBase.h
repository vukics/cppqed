// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_STRUCTURE_DYNAMICSBASE_H_INCLUDED
#define CPPQEDCORE_STRUCTURE_DYNAMICSBASE_H_INCLUDED

#include "DynamicsBaseFwd.h"

#include "ComplexExtensions.h"

#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>

#include <list>
#include <tuple>



namespace structure {


/// Provides services for dealing with frequency-like parameters, both real and complex, for all elements, frees and interactions alike, which are hence all derived from this class. 
/** 
 * Such parameters need a special treatment because every such parameter from all the subsystems (either frees or interactions) of the given physical system has to be considered
 * as a possible largest frequency of the whole system. This largest frequency is needed for determining the initial time-step of the ODE routine.
 * 
 * On the other hand, these parameters (together with all other, non-frequency-like parameters) have to be communicated towards the user when the framework
 * summarizes the parameters of a given run. Therefore, for each such parameter, the class stores not only the value, but also the name of the parameter,
 * plus another real number which multiplies the value of the named frequency, to give the actual frequency as appearing in the ODE.
 * 
 * An example will make this clear. Consider the Hamiltonian of a free mode: \f[\omega a^\dagger a.\f]
 * In this case, the parameter supplied by the user is \f$\omega\f$, but the largest frequency appearing in the ODE is actually this times the dimension of the system
 * (the cutoff of the Fock space). Hence, the tuple stored by the class for this particular frequency-like parameter will be something like:
 * ~~~
 * ("omega",omega,cutoff)
 * ~~~
 * 
 * In the case of a pumping term of the Hamiltonian: \f[\eta\lp a^\dagger+a\rp,\f] the multiplier will be different:
 * ~~~
 * "eta",eta,sqrt(cutoff)
 * ~~~
 * 
 * The class also stores an \refStdCppConstruct{std::ostringstream,sstream/ostringstream/} object, on which the constructor of any client (derived element class) can write its parameters,
 * and these will in turn be displayed when QuantumSystem::displayParameters is called for the system. Cf. \ref structurebundleguide "the structure-bundle guide".
 * 
 */
class DynamicsBase : private boost::noncopyable
{
public:
  typedef boost::shared_ptr<const DynamicsBase> Ptr;

  typedef std::tuple<std::string,double,double> RF; ///< name-value-multiplier tuple for a real frequency-like parameter
  typedef std::tuple<std::string,dcomp ,double> CF; ///< same for complex
  
  typedef std::list<RF>    RealFreqs; ///< list of real frequency-like parameters
  typedef std::list<CF> ComplexFreqs; ///< same for complex
  
  typedef std::initializer_list<RF>    RealFreqsInitializer;
  typedef std::initializer_list<CF> ComplexFreqsInitializer;

  static const    RealFreqs emptyRF;
  static const ComplexFreqs emptyCF;
  
  explicit DynamicsBase(const RealFreqs& =emptyRF, const ComplexFreqs& =emptyCF); ///< Straightforward constructor
  
  explicit DynamicsBase(const ComplexFreqs& complexFreqs) : DynamicsBase(emptyRF,complexFreqs) {}
  explicit DynamicsBase(RealFreqsInitializer rf, ComplexFreqsInitializer cf={}) : DynamicsBase(RealFreqs(rf),ComplexFreqs(cf)) {} ///< Constructor with initializer lists
  explicit DynamicsBase(ComplexFreqsInitializer cf) : DynamicsBase({},cf) {}
  explicit DynamicsBase(RF rf, CF cf=CF()) : DynamicsBase(RealFreqsInitializer{rf}, cf==CF() ? ComplexFreqsInitializer{} : ComplexFreqsInitializer{cf}) {}
  explicit DynamicsBase(CF cf) : DynamicsBase(ComplexFreqsInitializer{cf}) {}
           DynamicsBase(RealFreqsInitializer rf, CF cf) : DynamicsBase(rf,{cf}) {}
           DynamicsBase(RF rf, ComplexFreqsInitializer cf) : DynamicsBase({rf},cf) {}
  
  double highestFrequency() const; ///< Calculates the fastest timescale of the system from the frequencies stored in the lists

  std::ostream& displayParameters(std::ostream&) const; ///< Displays the content of the stored ostringstream followed by a call to #displayMoreParameters

  virtual ~DynamicsBase() {}

protected:
  std::ostringstream& getParsStream() {return paramsStream_;} ///< The stored `std::ostringstream` object, for constructors of clients to write parameters on
  
  virtual std::ostream& displayMoreParameters(std::ostream&) const; ///< In its default implementation, displayes the frequency-like parameters of the system in a nice format together with their names

  /// \name Getters
  //@{
     RealFreqs&    getRealFreqs() {return    realFreqs_;}
  ComplexFreqs& getComplexFreqs() {return complexFreqs_;}
  //@}
private:
  RealFreqs       realFreqs_;
  ComplexFreqs complexFreqs_;

  std::ostringstream paramsStream_;

};


} // structure



#endif // CPPQEDCORE_STRUCTURE_DYNAMICSBASE_H_INCLUDED
