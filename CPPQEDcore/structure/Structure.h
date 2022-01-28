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


template<typename T> auto castEx(std::shared_ptr<const T> qs) {return std::dynamic_pointer_cast<const Exact<T::N_RANK>>(qs);}
template<typename T> auto castHa(std::shared_ptr<const T> qs) {return std::dynamic_pointer_cast<const Hamiltonian<T::N_RANK>>(qs);}
template<typename T> auto castLi(std::shared_ptr<const T> qs) {return std::dynamic_pointer_cast<const Liouvillean<T::N_RANK>>(qs);}
template<typename T> auto castAv(std::shared_ptr<const T> qs) {return std::dynamic_pointer_cast<const Averaged<T::N_RANK>>(qs);}


/// Wrappers for Exact, Hamiltonian, Liouvillean, and Averaged member functions
/**
 * The aim is to determine whether the passed DynamicsBase or QuantumSystem object derives also from Exact, Hamiltonian, Liouvillean, or Averaged
 *
 * If the answer is
 * - positive, it forwards the member functions of the given class with the given object dynamic-cast to the necessary type.
 * - negative, it forwards the member functions of the given class as no-ops.
 *
 * \tparamRANK
 */
//@{


template<typename T>
std::ostream& streamCharacteristics(std::shared_ptr<const T> qs, std::ostream& os)
{
  return os<<"System characteristics: "
      <<(castEx(qs) ? "Interaction picture, "   : "")
      <<(castHa(qs) ? "Hamiltonian evolution, " : "")
      <<(castLi(qs) ? "Liouvillean evolution, " : "")
      <<(castAv(qs) ? "calculates Averages."    : "");
}


/// \name Forwarded members from Exact
//@{
template<typename T>
bool applicableInMaster(std::shared_ptr<const T> qs)
{
  if (const auto ex_=castEx(qs))
    return ex_->applicableInMaster();
  else return true;
}

template<typename T>
void actWithU(std::shared_ptr<const T> qs, double t, StateVectorLow<T::N_RANK>& psi, double t0)
{
  if (const auto ex_=castEx(qs)) ex_->actWithU(t,psi,t0);
}
//@}


/// \name Forwarded member from Hamiltonian
//@{
template<typename T>
void addContribution(std::shared_ptr<const T> qs, double t, const StateVectorLow<T::N_RANK>& psi, StateVectorLow<T::N_RANK>& dpsidt, double t0)
{
  if (const auto ha_=castHa(qs)) ha_->addContribution(t,psi,dpsidt,t0);
}
//@}


/// \name Forwarded member from Liouvillean
//@{
template<typename T>
void actWithJ(std::shared_ptr<const T> qs, double t, StateVectorLow<T::N_RANK>& psi, size_t lindbladNo)
{
  if (const auto li_=castLi(qs)) li_->actWithJ(t,psi,lindbladNo);
}

template<typename T>
void actWithSuperoperator(std::shared_ptr<const T> qs, double t, const DensityOperatorLow<T::N_RANK>& rho, DensityOperatorLow<T::N_RANK>& drhodt, size_t lindbladNo)
{
  if (const auto li_=castLi(qs)) li_->actWithSuperoperator(t,rho,drhodt,lindbladNo);
}
//@}

  
/// \name Forwarded members from Averaged
//@{
template<typename T>
void process(std::shared_ptr<const T> qs, Averages& averages)
{
  if (const auto av_=castAv(qs)) av_->process(averages);
}


/// If the first argument is a valid pointer, it calles Averaged::average, Averaged::process, and Averaged::stream in succession; otherwise a no-op. \related QuantumSystemWrapper
template<typename T>
std::tuple<std::ostream&,Averages>
stream(std::shared_ptr<const T> qs, double t, const quantumdata::LazyDensityOperator<T::N_RANK>& matrix, std::ostream& os, int precision)
{
  Averages averages;
  if (const auto av=castAv(qs)) {
    averages.reference(av->average(t,matrix));
    av->process(averages);
    av->stream(averages,os,precision);
  }
  return {os,averages};
}

//@}


/// \name Forwarded members from LiouvilleanAveragedCommon
//@{

template<LiouvilleanAveragedTag LA, typename T>
std::shared_ptr<const LiouvilleanAveragedCommonRanked<T::N_RANK>> dispatchLiouvilleanAverage(std::shared_ptr<const T> qs)
{
  if constexpr (LA==LA_Li) return castLi(qs);
  else return castAv(qs);
}


template<LiouvilleanAveragedTag LA, typename T>
size_t nAvr(std::shared_ptr<const T> qs)
{
  const auto ptr=dispatchLiouvilleanAverage<LA>(qs);
  return ptr ? ptr->nAvr() : 0;
}

template<LiouvilleanAveragedTag LA, typename T>
std::ostream& streamKey(std::shared_ptr<const T> qs, std::ostream& os, size_t& i)
{
  if (const auto ptr=dispatchLiouvilleanAverage<LA>(qs)) ptr->streamKey(os,i);
  return os;
}


template<LiouvilleanAveragedTag LA, typename T>
const Averages average(std::shared_ptr<const T> qs, double t, const quantumdata::LazyDensityOperator<T::N_RANK>& matrix)
{
  if (const auto ptr=dispatchLiouvilleanAverage<LA>(qs))
    return ptr->average(t,matrix);
  else return Averages{};
}
//@}


//@}


} // structure



#endif // CPPQEDCORE_STRUCTURE_STRUCTURE_H_INCLUDED
