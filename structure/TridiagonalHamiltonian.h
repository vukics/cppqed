// -*- C++ -*-
/// \briefFileDefault
#ifndef STRUCTURE_TRIDIAGONALHAMILTONIAN_H_INCLUDED
#define STRUCTURE_TRIDIAGONALHAMILTONIAN_H_INCLUDED

#include "TridiagonalHamiltonianFwd.h"

#include "Hamiltonian.h"

#include "Tridiagonal.h"

#include <list>


namespace structure {


namespace details {


/// TridiagonalHamiltonian base
template<int RANK, bool IS_TIME_DEPENDENT>
class TDH_Base 
  : public Hamiltonian<RANK,
                       mpl::if_c<IS_TIME_DEPENDENT,
                                 mpl::integral_c<TimeDependence,ONE_TIME>,
                                 mpl::integral_c<TimeDependence, NO_TIME>
                                 >::type::value
                       >
{
protected:
  typedef quantumoperator::Tridiagonal<RANK> Tridiagonal ;
  typedef std::list<Tridiagonal>             Tridiagonals;

  typedef typename Hamiltonian<RANK>::StateVectorLow StateVectorLow; ///< Should be the same for Hamiltonian classes of any degree of TimeDependence.

  TDH_Base(const Tridiagonals& hOverIs) : hOverIs_(hOverIs) {}

  const Tridiagonals& getH_OverIs() const {return hOverIs_;}
        Tridiagonals& getH_OverIs()       {return const_cast<Tridiagonals&>(static_cast<const TDH_Base*>(this)->getH_OverIs());}

private:
  /// \name Virtuals inherited from Hamiltonian
  //@{
    /**
     * The first implements Hamiltonian<RANK,ONE_TIME> and can be used only when `IS_TIME_DEPENDENT=true`, while the second implements Hamiltonian<RANK,NO_TIME>
     * and can be used only in the opposite case.
     * 
     * \note The technique of providing both functions here depends on that feature of C++ templates that only such parts of a template are instantiated as are actually invoked.
     * This could not work if Hamiltonian<RANK,NO_TIME> inherited from Hamiltonian<RANK,ONE_TIME> instead of Hamiltonian<RANK,TWO_TIME>, because in that case, both functions
     * would be implementations of virtual functions and hence (and the first one, in the case of `IS_TIME_DEPENDENT=false` would even silently re-implement the virtual function
     * inherited from Hamiltonian<RANK,ONE_TIME> through Hamiltonian<RANK,NO_TIME>) *always instantiated*!
     * 
     * \see the technique used @ FreeExact & ElementLiouvillean
     * 
     */
  void addContribution_v(double, const StateVectorLow&, StateVectorLow&) const;
  void addContribution_v(        const StateVectorLow&, StateVectorLow&) const;
  //@}
  mutable Tridiagonals hOverIs_;

};


} // details



/// Implements the action of a Hamiltonian whose matrix consists of a sum of \link quantumoperator::Tridiagonal tridiagonal matrices\endlink
/**
 * \f[H_\text{tridiagonals}(t)=H_0(t)+H_1(t)+H_2(t)+\dots\f] with the \f$H_i(t)\f$s being all described by quantumoperator::Tridiagonal `<RANK>` objects.
 * 
 * Such a class can be constructed with either a list of quantumoperator::Tridiagonal `<RANK>` objects, or only one such object when the above sum consists of only one term.
 * 
 * \tparamRANK
 * \tparam IS_TIME_DEPENDENT governs time-dependence & the composition of the class @ compile time
 * 
 * Implements Hamiltonian<RANK,ONE_TIME> when `IS_TIME_DEPENDENT=true` **OR** Hamiltonian<RANK,NO_TIME>  when `IS_TIME_DEPENDENT=false`
 * 
 * \note The present architecture of quantumoperator::Tridiagonal does not allow to cover the case #TWO_TIME.
 *  
 */
template<int RANK, bool IS_TIME_DEPENDENT>
class TridiagonalHamiltonian : public details::TDH_Base<RANK,IS_TIME_DEPENDENT>
{
private:
  typedef details::TDH_Base<RANK,IS_TIME_DEPENDENT> Base;

public:
  typedef typename Base::Tridiagonal  Tridiagonal ;
  typedef typename Base::Tridiagonals Tridiagonals;

  TridiagonalHamiltonian(const Tridiagonals& hOverIs) : Base(hOverIs) {} ///< Generic constructor
  TridiagonalHamiltonian(const Tridiagonal & hOverI ) : Base(Tridiagonals(1,hOverI)) {} ///< \overload

};

#undef BASE_class

} // structure


#endif // STRUCTURE_TRIDIAGONALHAMILTONIAN_H_INCLUDED
