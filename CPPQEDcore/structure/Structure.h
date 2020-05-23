// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines structur::QuantumSystemWrapper}
#ifndef CPPQEDCORE_STRUCTURE_STRUCTURE_H_INCLUDED
#define CPPQEDCORE_STRUCTURE_STRUCTURE_H_INCLUDED

#include "QuantumSystem.h"
#include "DynamicsBase.h"

#include "Exact.h"
#include "Hamiltonian.h"
#include "Liouvillean.h"
#include "Averaged.h"

/// Comprises modules for describing quantum systems.
/** 
 * Among them the most important is QuantumSystem, which is an abstract interface every system has to provide to be usable with quantum trajectories like quantumtrajectory::MCWF_Trajectory 
 * or quantumtrajectory::Master. This is why all the elementary and composite systems are more or less directly derived from QuantumSystem.
 * 
 * Much of the design here depends on the requirements of a step of the Monte-Carlo wave function method, as described in \ref mcwftrajectory, so the reader is asked to have a look at there, too.
 * 
 * Most of the classes in this namespace belong to a single hierarchy, sketched below. This diagram is however approximate, because the picture in reality is somewhat complicated 
 * by the heavy use of templates, partial specializations, and conditional inheritance. The broken lines signifying that the inheritance is not direct, due to some classes in between,
 * which can be considered implementation details.
 * 
 * \image html structure.png
 * 
 * We have also indicated how classes representing elementary free subsystems (Mode) and interactions (JaynesCummings), and those representing composite systems (BinarySystem and Composite)
 * fit into the hierarchy.
 *
 * \note This figure is slightly outdated, but still represents the basic idea
 * 
 * These modules in the hierarchy provide a lot of services for implementing new elements in the framework. For examples on how to optimally use these services, 
 * \ref structurebundleguide "the structure-bundle guide".
 * 
 */
namespace structure {


typedef blitz::TinyVector<bool,3> SystemCharacteristics;


using std::dynamic_pointer_cast;


/// Dynamic cast to a shared pointer to Exact
template<int RANK, typename T>
inline auto
qse(std::shared_ptr<const T> t)
{return dynamic_pointer_cast<const Exact<RANK> >(t);}

/// Dynamic cast to a shared pointer to Hamiltonian
template<int RANK, typename T>
inline auto
qsh(std::shared_ptr<const T> t)
{return dynamic_pointer_cast<const Hamiltonian<RANK> >(t);}

/// Dynamic cast to a shared pointer to Liouvillean
template<int RANK, typename T>
inline auto
qsl(std::shared_ptr<const T> t)
{return dynamic_pointer_cast<const Liouvillean<RANK> >(t);}

/// Dynamic cast to a shared pointer to Averaged
template<int RANK, typename T>
inline auto
qsa(std::shared_ptr<const T> t)
{return dynamic_pointer_cast<const Averaged<RANK> >(t);}


// Some functions that are used in contexts other than QuantumSystemWrapper are factored out:


/// If the first argument is a valid pointer, it calles Averaged::average, Averaged::process, and Averaged::display in succession; otherwise a no-op. \related QuantumSystemWrapper
template<int RANK>
std::ostream& display(std::shared_ptr<const Averaged<RANK> > av,
                      double t,
                      const quantumdata::LazyDensityOperator<RANK>& matrix,
                      std::ostream& os,
                      int precision)
{
  if (av) {
    typename Averaged<RANK>::Averages averages(av->average(t,matrix));
    av->process(averages);
    av->display(averages,os,precision);
  }
  return os;
}


/// If the first argument is a valid pointer, it calles LiouvilleanAveragedCommon Averaged::average; otherwise a no-op (returning an empty array) \related QuantumSystemWrapper
template<int RANK>
const LiouvilleanAveragedCommon::DArray1D average(typename LiouvilleanAveragedCommonRanked<RANK>::Ptr ptr, double t, const quantumdata::LazyDensityOperator<RANK>& matrix)
{
  return ptr ? ptr->average(t,matrix) : LiouvilleanAveragedCommon::DArray1D();
}


/// A wrapper for Exact, Hamiltonian, Liouvillean, and Averaged
/**
 * It’s aim is to determine whether the passed DynamicsBase or QuantumSystem object derives also from Exact, Hamiltonian, Liouvillean, and Averaged
 *
 * If the answer is
 * - positive, it forwards the member functions of the given class with the given object dynamic-cast to the necessary type.
 * - negative, it forwards the member functions of the given class as no-ops.
 *
 * \tparamRANK
 * \tparam IS_CONST governs const-ness (non-const necessary for assignment)
 *
 */
template<int RANK, bool IS_CONST> 
class QuantumSystemWrapper
{
public:
  static const int N_RANK=RANK;
  
  /// \name Wrapped types
  //@{
  typedef QuantumSystem<RANK> QS;
  typedef Exact        <RANK> Ex;
  typedef Hamiltonian  <RANK> Ha;
  typedef Liouvillean  <RANK> Li;
  typedef Averaged     <RANK> Av;
  //@}

  /// \name Pointers to wrapped types
  //@{
  typedef typename QS::Ptr QuantumSystemPtr;
  typedef typename Ex::Ptr         ExactPtr;
  typedef typename Ha::Ptr   HamiltonianPtr;
  typedef typename Li::Ptr   LiouvilleanPtr;
  typedef typename Av::Ptr      AveragedPtr;
  //@}

  /// \name Necessary types from wrapped types for definition of member-function signatures
  //@{
  typedef typename Ex::StateVectorLow StateVectorLow;

  typedef typename Li::Rates               Rates              ;
  typedef typename Li::LazyDensityOperator LazyDensityOperator;

  typedef typename Li::DensityOperatorLow DensityOperatorLow;
  
  typedef typename Av::Averages Averages;
  //@}

  /// \name Constructors
  //@{
  explicit QuantumSystemWrapper(DynamicsBase::Ptr qs)
    : qs_(dynamic_pointer_cast<const QuantumSystem<RANK> >(qs)),
      ex_(qse<RANK>(qs)),
      ha_(qsh<RANK>(qs)),
      li_(qsl<RANK>(qs)),
      av_(qsa<RANK>(qs))
  {} ///< Constructor from DynamicsBase

  explicit QuantumSystemWrapper(QuantumSystemPtr qs, bool isNoisy)
    : qs_(qs),
      ex_(qse<RANK>(qs)),
      ha_(qsh<RANK>(qs)),
      li_(isNoisy ? qsl<RANK>(qs) : LiouvilleanPtr()),
      av_(qsa<RANK>(qs))
  {} ///< Constructor from QuantumSystem
  //@}

  /// \name Getters
  //@{
  QuantumSystemPtr getQS() const {return qs_;}
  ExactPtr         getEx() const {return ex_;} 
  HamiltonianPtr   getHa() const {return ha_;}
  LiouvilleanPtr   getLi() const {return li_;} 
  AveragedPtr      getAv() const {return av_;}  
  //@}

private:
  typedef typename LiouvilleanAveragedCommonRanked<RANK>::Ptr L_or_A_Ptr;

public:
  /**
   * \name Dispatcher between Liouvillean and Averaged
   * \internal We use overload instead of template specialization, which is only possible in namespace scope
   *  @{
   */
  const L_or_A_Ptr getLA(LA_Li_tagType) const {return li_;}
  const L_or_A_Ptr getLA(LA_Av_tagType) const {return av_;}
  //@}

  /// Displays the dynamical characteristics of the system
  std::ostream& displayCharacteristics(std::ostream& os) const {return os<<"System characteristics: "<<(ex_ ? "Interaction picture, "   : "")<<(ha_ ? "Hamiltonian evolution, " : "")<<(li_ ? "Liouvillean evolution, " : "")<<(av_ ? "calculates Averages."    : "");}


  /// \name Forwarded members from Exact
  //@{
  bool applicableInMaster() const {return ex_ ? ex_->applicableInMaster() : true;}

  void actWithU(double t, StateVectorLow& psi, double t0) const {if (ex_) ex_->actWithU(t,psi,t0);}
  //@}


  /// \name Forwarded member from Hamiltonian
  //@{
  void addContribution(double t, const StateVectorLow& psi, StateVectorLow& dpsidt, double t0) const {if (ha_) ha_->addContribution(t,psi,dpsidt,t0);}
  //@}


  /// \name Forwarded member from Liouvillean
  //@{
  void actWithJ(double t, StateVectorLow& psi, size_t lindbladNo) const {if (li_) li_->actWithJ(t,psi,lindbladNo);}
  void actWithSuperoperator(double t, const DensityOperatorLow& rho, DensityOperatorLow& drhodt, size_t lindbladNo) const {if (li_) li_->actWithSuperoperator(t,rho,drhodt,lindbladNo);}
  //@}

  
  /// \name Forwarded members from Averaged
  //@{
  void process(Averages& averages) const {if (av_) av_->process(averages);}

  std::ostream& display(double t, const LazyDensityOperator& matrix, std::ostream& os, int precision) const {return structure::display(av_,t,matrix,os,precision);}
  //@}


  /// \name Forwarded members from LiouvilleanAveragedCommon
  //@{
  template<LiouvilleanAveragedTag LA>
  size_t nAvr() const {const auto ptr=getLA(LiouvilleanAveragedTag_<LA>()); return ptr ? ptr->nAvr() : 0;}

  template<LiouvilleanAveragedTag LA>
  std::ostream& displayKey(std::ostream& os, size_t& i) const {if (const auto ptr=getLA(LiouvilleanAveragedTag_<LA>())) ptr->displayKey(os,i); return os;}

  template<LiouvilleanAveragedTag LA>
  const Averages average(double t, const LazyDensityOperator& matrix) const {return structure::average(getLA(LiouvilleanAveragedTag_<LA>()),t,matrix);}
  //@}

protected:
  QuantumSystemWrapper() : qs_(), ex_(), ha_(), li_(), av_() {}

private:
  typename tmptools::ConditionalAddConst<QuantumSystemPtr,IS_CONST>::type qs_;
  typename tmptools::ConditionalAddConst<ExactPtr        ,IS_CONST>::type ex_; 
  typename tmptools::ConditionalAddConst<HamiltonianPtr  ,IS_CONST>::type ha_;
  typename tmptools::ConditionalAddConst<LiouvilleanPtr  ,IS_CONST>::type li_; 
  typename tmptools::ConditionalAddConst<AveragedPtr     ,IS_CONST>::type av_;
  
};



} // structure



#endif // CPPQEDCORE_STRUCTURE_STRUCTURE_H_INCLUDED
