// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
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


// Some functions that are used in contexts other than QuantumSystemWrapper are factored out:

/// If the first argument is a valid pointer, it calles Averaged::average, Averaged::process, and Averaged::stream in succession; otherwise a no-op. \related QuantumSystemWrapper
template<int RANK>
std::tuple<std::ostream&,Averages>
stream(AveragedPtr<RANK> av,
        double t,
        const quantumdata::LazyDensityOperator<RANK>& matrix,
        std::ostream& os,
        int precision)
{
  Averages averages;
  if (av) {
    averages.reference(av->average(t,matrix));
    av->process(averages);
    av->stream(averages,os,precision);
  }
  return {os,averages};
}


/// If the first argument is a valid pointer, it calles LiouvilleanAveragedCommon Averaged::average; otherwise a no-op (returning an empty array) \related QuantumSystemWrapper
template<int RANK>
auto average(std::shared_ptr<const LiouvilleanAveragedCommonRanked<RANK>> ptr, double t, const quantumdata::LazyDensityOperator<RANK>& matrix)
{
  return ptr ? ptr->average(t,matrix) : Averages{};
}


/// A wrapper for Exact, Hamiltonian, Liouvillean, and Averaged
/**
 * Its aim is to determine whether the passed DynamicsBase or QuantumSystem object derives also from Exact, Hamiltonian, Liouvillean, and Averaged
 *
 * If the answer is
 * - positive, it forwards the member functions of the given class with the given object dynamic-cast to the necessary type.
 * - negative, it forwards the member functions of the given class as no-ops.
 *
 * \tparamRANK
 */
template<int RANK> 
class QuantumSystemWrapper
{
public:
  static const int N_RANK=RANK;
  
  /// \name Constructors
  //@{
  explicit QuantumSystemWrapper(DynamicsBasePtr qs)
    : qs_(std::dynamic_pointer_cast<const QuantumSystem<RANK>>(qs)),
      ex_(std::dynamic_pointer_cast<const Exact<RANK>>(qs)),
      ha_(std::dynamic_pointer_cast<const Hamiltonian<RANK>>(qs)),
      li_(std::dynamic_pointer_cast<const Liouvillean<RANK>>(qs)),
      av_(std::dynamic_pointer_cast<const Averaged<RANK>>(qs))
  {} ///< Constructor from DynamicsBase

  explicit QuantumSystemWrapper(QuantumSystemPtr<RANK> qs, bool isNoisy)
    : qs_(qs),
      ex_(std::dynamic_pointer_cast<const Exact<RANK>>(qs)),
      ha_(std::dynamic_pointer_cast<const Hamiltonian<RANK>>(qs)),
      li_(isNoisy ? std::dynamic_pointer_cast<const Liouvillean<RANK>>(qs) : LiouvilleanPtr<RANK>()),
      av_(std::dynamic_pointer_cast<const Averaged<RANK>>(qs))
  {} ///< Constructor from QuantumSystem
  //@}

  /// \name Getters
  //@{
  auto getQS() const {return qs_;}
  
  auto getEx() const {return ex_;} 
  auto getHa() const {return ha_;}
  auto getLi() const {return li_;} 
  auto getAv() const {return av_;}  
  //@}

public:
  /**
   * \name Dispatcher between Liouvillean and Averaged
   * \internal We use overload instead of template specialization, which is only possible in namespace scope
   *  @{
   */
  std::shared_ptr<const LiouvilleanAveragedCommonRanked<RANK>> getLA(LA_Li_tagType) const {return li_;}
  std::shared_ptr<const LiouvilleanAveragedCommonRanked<RANK>> getLA(LA_Av_tagType) const {return av_;}
  //@}

  /// Streams the dynamical characteristics of the system
  std::ostream& streamCharacteristics(std::ostream& os) const {return os<<"System characteristics: "<<(ex_ ? "Interaction picture, "   : "")<<(ha_ ? "Hamiltonian evolution, " : "")<<(li_ ? "Liouvillean evolution, " : "")<<(av_ ? "calculates Averages."    : "");}


  /// \name Forwarded members from Exact
  //@{
  bool applicableInMaster() const {return ex_ ? ex_->applicableInMaster() : true;}

  void actWithU(double t, StateVectorLow<RANK>& psi, double t0) const {if (ex_) ex_->actWithU(t,psi,t0);}
  //@}


  /// \name Forwarded member from Hamiltonian
  //@{
  void addContribution(double t, const StateVectorLow<RANK>& psi, StateVectorLow<RANK>& dpsidt, double t0) const {if (ha_) ha_->addContribution(t,psi,dpsidt,t0);}
  //@}


  /// \name Forwarded member from Liouvillean
  //@{
  void actWithJ(double t, StateVectorLow<RANK>& psi, size_t lindbladNo) const {if (li_) li_->actWithJ(t,psi,lindbladNo);}
  void actWithSuperoperator(double t, const DensityOperatorLow<RANK>& rho, DensityOperatorLow<RANK>& drhodt, size_t lindbladNo) const {if (li_) li_->actWithSuperoperator(t,rho,drhodt,lindbladNo);}
  //@}

  
  /// \name Forwarded members from Averaged
  //@{
  void process(Averages& averages) const {if (av_) av_->process(averages);}

  std::tuple<std::ostream&,Averages> stream(double t, const quantumdata::LazyDensityOperator<RANK>& matrix, std::ostream& os, int precision) const
  {
    return structure::stream(av_,t,matrix,os,precision);
  }
  //@}


  /// \name Forwarded members from LiouvilleanAveragedCommon
  //@{
  template<LiouvilleanAveragedTag LA>
  size_t nAvr() const {const auto ptr=getLA(LiouvilleanAveragedTag_<LA>()); return ptr ? ptr->nAvr() : 0;}

  template<LiouvilleanAveragedTag LA>
  std::ostream& streamKey(std::ostream& os, size_t& i) const {if (const auto ptr=getLA(LiouvilleanAveragedTag_<LA>())) ptr->streamKey(os,i); return os;}

  template<LiouvilleanAveragedTag LA>
  const Averages average(double t, const quantumdata::LazyDensityOperator<RANK>& matrix) const {return structure::average(getLA(LiouvilleanAveragedTag_<LA>()),t,matrix);}
  //@}

protected:
  QuantumSystemWrapper() : qs_(), ex_(), ha_(), li_(), av_() {}

private:
  QuantumSystemPtr<RANK> qs_;
  
  ExactPtr<RANK> ex_;

  HamiltonianPtr<RANK> ha_;
  
  LiouvilleanPtr<RANK> li_;
  
  AveragedPtr<RANK> av_;
  
};



} // structure



#endif // CPPQEDCORE_STRUCTURE_STRUCTURE_H_INCLUDED
