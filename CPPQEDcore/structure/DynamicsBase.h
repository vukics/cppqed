// -*- C++ -*-
/// \briefFileDefault
#ifndef STRUCTURE_DYNAMICSBASE_H_INCLUDED
#define STRUCTURE_DYNAMICSBASE_H_INCLUDED

#include "DynamicsBaseFwd.h"

#include "ComplexExtensions.h"

#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>
#include <boost/tuple/tuple.hpp>

#include <boost/assign/list_of.hpp>

#include <list>


/// Macro facilitating the creation of lists of the name-value-multiplier tuples
/** Cf. \refBoost{the <tt>tuple_list_of</tt> tool from  Boost.Assign,assign/doc/index.html#tuple_list_of} */
#define FREQS boost::assign::tuple_list_of


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
 * The class also stores an \refStdCppConstruct{std::stringstream,sstream/stringstream/} object, on which the constructor of any client (derived element class) can write its parameters,
 * and these will in turn be displayed when QuantumSystem::displayParameters is called for the system. Cf. \ref structurebundleguide "the structure-bundle guide".
 * 
 */
class DynamicsBase : private boost::noncopyable
{
public:
  typedef boost::shared_ptr<const DynamicsBase> Ptr;

  typedef std::list<boost::tuple<std::string,double,double> >    RealFreqs; ///< Stores the name-value-multiplier \refBoost{tuples,tuple/doc/tuple_users_guide.html} for real frequency-like parameters
  typedef std::list<boost::tuple<std::string,dcomp ,double> > ComplexFreqs; ///< same for complex

  explicit DynamicsBase(const    RealFreqs& =RealFreqs   (), 
                        const ComplexFreqs& =ComplexFreqs()); ///< Straightforward constructor

  double highestFrequency() const; ///< Calculates the fastest timescale of the system from the frequencies stored in the lists

  std::ostream& displayParameters(std::ostream&) const; ///< Displays the content of the stored stringstream followed by a call to #displayMoreParameters

  virtual ~DynamicsBase() {}

protected:
  std::stringstream& getParsStream() {return paramsStream_;} ///< The stored `std::stringstream` object, for constructors of clients to write parameters on
  
  virtual std::ostream& displayMoreParameters(std::ostream&) const; ///< In its default implementation, displayes the frequency-like parameters of the system in a nice format together with their names

  /// \name Getters
  //@{
     RealFreqs&    getRealFreqs() {return    realFreqs_;}
  ComplexFreqs& getComplexFreqs() {return complexFreqs_;}
  //@}  
private:
  RealFreqs       realFreqs_;
  ComplexFreqs complexFreqs_;

  std::stringstream paramsStream_;

};


} // structure



#endif // STRUCTURE_DYNAMICSBASE_H_INCLUDED
