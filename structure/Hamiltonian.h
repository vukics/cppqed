// -*- C++ -*-
/// \briefFile{Defines the hierarchical partial specializations of structure::Hamiltonian}
#ifndef STRUCTURE_HAMILTONIAN_H_INCLUDED
#define STRUCTURE_HAMILTONIAN_H_INCLUDED

#include "HamiltonianFwd.h"

#include "Types.h"


namespace structure {


/// The first partial specialization of the general template Hamiltonian for the most general degree of \link TimeDependence time dependence\endlink (#TWO_TIME).
/** \tparamRANK */
template<int RANK>
class Hamiltonian<RANK,TWO_TIME> : private quantumdata::Types<RANK>
{
public:
  typedef boost::shared_ptr<const Hamiltonian> Ptr;

  typedef typename quantumdata::Types<RANK>::StateVectorLow StateVectorLow;

  /// Adds the Hamiltonian contribution \f$\frac{H(t)}i\ket\Psi\f$ of the given (sub)system to `dpsidt`
  /**
   * The assumption is that the time when the Schrödinger picture and interaction picture (if any) coincide is `tIntPic0`. There are two important points to note:
   * -# The contribution has to be *added* to `dpsidt` instead of `dpsidt` being *replaced*. 
   * This is because when the given system is embedded in a larger system, other (sub)systems may also contribute.
   * -# The function has to calculate the effect of \f$\frac{H(t)}i\f$ and not merely \f$H\f$, since it is the former which determines the derivative of the state vector.
   * 
   * This latter is so often missed, that we emphasize it again (although we know that it will still be missed from time to time):
   * \warning When implementing the Hamiltonian, not \f$H\f$ itself but \f$\frac Hi\f$ has to supplied!
   * 
   */
  void addContribution(double t, ///<[in] the time instant \f$t\f$ for #TWO_TIME dependence
                       const StateVectorLow& psi, ///<[in] the state vector \f$\ket\Psi\f$
                       StateVectorLow& dpsidt, ///<[in,out] the state vector to be contributed to by \f$\frac{H(t)}i\ket\Psi\f$
                       double tIntPic0 ///<[in] the time instant \f$t_0\f$ for #TWO_TIME dependence
                      ) const
                       {addContribution_v(t,psi,dpsidt,tIntPic0);}

  virtual ~Hamiltonian() {}

private:
  virtual void addContribution_v(double, const StateVectorLow&, StateVectorLow&, double) const = 0; 

};


/// The second partial specialization of the general template Hamiltonian for the middle degree of \link TimeDependence time dependence\endlink (#ONE_TIME).
template<int RANK>
class Hamiltonian<RANK,ONE_TIME> : public Hamiltonian<RANK>
{
public:
  typedef typename Hamiltonian<RANK>::StateVectorLow StateVectorLow;

private:
  /// The inherited virtual gets implemented by calling a newly defined virtual with only one time parameter
  void addContribution_v(double t, const StateVectorLow& psi, StateVectorLow& dpsidt, double tIntPic0) const
  {
    addContribution_v(t-tIntPic0,psi,dpsidt);
  }

  /// The newly defined virtual with only one time parameter
  virtual void addContribution_v(double, const StateVectorLow&, StateVectorLow&) const = 0;
  
};


/// The third partial specialization of the general template Hamiltonian for \link TimeDependence no time dependence\endlink (#NO_TIME).
template<int RANK>
class Hamiltonian<RANK,NO_TIME> : public Hamiltonian<RANK>
{
public:
  typedef typename Hamiltonian<RANK>::StateVectorLow StateVectorLow;

private:
  /// The inherited virtual gets implemented by calling a newly defined virtual with no time parameter
  void addContribution_v(double, const StateVectorLow& psi, StateVectorLow& dpsidt, double) const 
  {
    addContribution_v(psi,dpsidt);
  }

  /// The newly defined virtual with no time parameter
  virtual void addContribution_v(const StateVectorLow&, StateVectorLow&) const = 0;
  
};



} // structure

#endif // STRUCTURE_HAMILTONIAN_H_INCLUDED
