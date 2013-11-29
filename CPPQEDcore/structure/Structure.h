// -*- C++ -*-
#ifndef STRUCTURE_STRUCTURE_H_INCLUDED
#define STRUCTURE_STRUCTURE_H_INCLUDED

#include "StructureFwd.h"

#include "QuantumSystem.h"
#include "DynamicsBase.h"

#include "Exact.h"
#include "Hamiltonian.h"
#include "Liouvillean.h"
#include "Averaged.h"


// Note that, in each case, the most general class is taken.
// This is because Hamiltonian ONE_TIME & NO_TIME is derived from TWO_TIME, Liouvillean false is derived from true, and Averaged false is derived from true.
// At the same time, these most general cases are the default ones.


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
 * These modules in the hierarchy provide a lot of services for implementing new elements in the framework. For examples on how to optimally use these services, 
 * \ref structurebundleguide "the structure-bundle guide".
 * 
 */
namespace structure {


typedef blitz::TinyVector<bool,3> SystemCharacteristics;


using boost::dynamic_pointer_cast;

template<int RANK>
inline 
const typename Exact<RANK>::Ptr 
qse(boost::shared_ptr<const QuantumSystem<RANK> > quantumSystem)
{return dynamic_pointer_cast<const Exact<RANK> >(quantumSystem);}

template<int RANK>
inline 
const typename Hamiltonian<RANK>::Ptr 
qsh(boost::shared_ptr<const QuantumSystem<RANK> > quantumSystem)
{return dynamic_pointer_cast<const Hamiltonian<RANK> >(quantumSystem);}

template<int RANK>
inline 
const typename Liouvillean<RANK>::Ptr 
qsl(boost::shared_ptr<const QuantumSystem<RANK> > quantumSystem)
{return dynamic_pointer_cast<const Liouvillean<RANK> >(quantumSystem);}

template<int RANK>
inline 
const typename Averaged<RANK>::Ptr 
qsa(boost::shared_ptr<const QuantumSystem<RANK> > quantumSystem)
{return dynamic_pointer_cast<const Averaged<RANK> >(quantumSystem);}



template<int RANK>
inline 
const typename Exact<RANK>::Ptr 
qse(DynamicsBase::Ptr base)
{return dynamic_pointer_cast<const Exact<RANK> >(base);}

template<int RANK>
inline 
const typename Hamiltonian<RANK>::Ptr 
qsh(DynamicsBase::Ptr base)
{return dynamic_pointer_cast<const Hamiltonian<RANK> >(base);}

template<int RANK>
inline 
const typename Liouvillean<RANK>::Ptr 
qsl(DynamicsBase::Ptr base)
{return dynamic_pointer_cast<const Liouvillean<RANK> >(base);}

template<int RANK>
inline 
const typename Averaged<RANK>::Ptr 
qsa(DynamicsBase::Ptr base)
{return dynamic_pointer_cast<const Averaged<RANK> >(base);}


// Some functions that are used in contexts other than QuantumSystemWrapper are factored out:

template<int RANK>
std::ostream& display(boost::shared_ptr<const Averaged<RANK> >, double, const quantumdata::LazyDensityOperator<RANK>&, std::ostream&, int);



template<int RANK>
const LiouvilleanAveragedCommon::DArray1D average(typename LiouvilleanAveragedCommonRanked<RANK>::Ptr, double, const quantumdata::LazyDensityOperator<RANK>&);



template<int RANK, bool IS_CONST> 
class QuantumSystemWrapper
{
public:
  static const int N_RANK=RANK;
  
  typedef QuantumSystem<RANK> QS;
  typedef Exact        <RANK> Ex;
  typedef Hamiltonian  <RANK> Ha;
  typedef Liouvillean  <RANK> Li;
  typedef Averaged     <RANK> Av;

  typedef typename QS::Ptr QuantumSystemPtr;
  typedef typename Ex::Ptr         ExactPtr;
  typedef typename Ha::Ptr   HamiltonianPtr;
  typedef typename Li::Ptr   LiouvilleanPtr;
  typedef typename Av::Ptr      AveragedPtr;

  typedef typename Ex::StateVectorLow StateVectorLow;

  typedef typename Li::Rates               Rates              ;
  typedef typename Li::LazyDensityOperator LazyDensityOperator;

  typedef typename Av::Averages Averages;

  explicit QuantumSystemWrapper(DynamicsBase::Ptr qs)
    : qs_(dynamic_pointer_cast<const QuantumSystem<RANK> >(qs)),
      ex_(qse<RANK>(qs)),
      ha_(qsh<RANK>(qs)),
      li_(qsl<RANK>(qs)),
      av_(qsa<RANK>(qs))
  {}

  explicit QuantumSystemWrapper(QuantumSystemPtr qs, bool isNoisy=true)
    : qs_(qs),
      ex_(qse<RANK>(qs)),
      ha_(qsh<RANK>(qs)),
      li_(isNoisy ? qsl<RANK>(qs) : LiouvilleanPtr()),
      av_(qsa<RANK>(qs))
  {}

  const QuantumSystemPtr getQS() const {return qs_;}
  const ExactPtr         getEx() const {return ex_;} 
  const HamiltonianPtr   getHa() const {return ha_;}
  const LiouvilleanPtr   getLi() const {return li_;} 
  const AveragedPtr      getAv() const {return av_;}  

  QuantumSystemPtr getQS() {return qs_;}
  ExactPtr         getEx() {return ex_;} 
  HamiltonianPtr   getHa() {return ha_;}
  LiouvilleanPtr   getLi() {return li_;} 
  AveragedPtr      getAv() {return av_;}  

private:
  typedef typename LiouvilleanAveragedCommonRanked<RANK>::Ptr L_or_A_Ptr;

public:  
  // overload instead of template specialization, which is only possible in namespace scope
  const L_or_A_Ptr getLA(LA_Li_tagType) const {return li_;}
  const L_or_A_Ptr getLA(LA_Av_tagType) const {return av_;}
  

  std::ostream& displayCharacteristics(std::ostream& os) const {return os<<"# System characteristics: "<<(ex_ ? "Interaction picture, "   : "")<<(ha_ ? "Hamiltonian evolution, " : "")<<(li_ ? "Liouvillean evolution, " : "")<<(av_ ? "calculates Averages."    : "");}

  
  // Exact
  
  bool isUnitary() const {return ex_ ? ex_->isUnitary() : true;}

  void actWithU(double t, StateVectorLow& psi, double t0) const {if (ex_) ex_->actWithU(t,psi,t0);}


  // Hamiltonian
  
  void addContribution(double t, const StateVectorLow& psi, StateVectorLow& dpsidt, double tIntPic0) const {if (ha_) ha_->addContribution(t,psi,dpsidt,tIntPic0);}


  // Liouvillean
  
  void actWithJ(double t, StateVectorLow& psi, size_t lindbladNo) const {if (li_) li_->actWithJ(t,psi,lindbladNo);}

  
  // Averaged
  
  void process(Averages& averages) const {if (av_) av_->process(averages);}

  std::ostream& display(double t, const LazyDensityOperator& matrix, std::ostream& os, int precision) const {return structure::display(av_,t,matrix,os,precision);}


  // LiouvilleanAveragedCommon

  template<LiouvilleanAveragedTag LA>
  size_t nAvr() const {const L_or_A_Ptr ptr=getLA(LiouvilleanAveragedTag_<LA>()); return ptr ? ptr->nAvr() : 0;}

  template<LiouvilleanAveragedTag LA>
  std::ostream& displayKey(std::ostream& os, size_t& i) const {if (const L_or_A_Ptr ptr=getLA(LiouvilleanAveragedTag_<LA>())) ptr->displayKey(os,i); return os;}

  template<LiouvilleanAveragedTag LA>
  const Averages average(double t, const LazyDensityOperator& matrix) const {return structure::average(getLA(LiouvilleanAveragedTag_<LA>()),t,matrix);}
  
protected:
  QuantumSystemWrapper() : qs_(), ex_(), ha_(), li_(), av_() {}

private:
  typename tmptools::ConditionalAddConst<QuantumSystemPtr,IS_CONST>::type qs_;
  typename tmptools::ConditionalAddConst<ExactPtr        ,IS_CONST>::type ex_; 
  typename tmptools::ConditionalAddConst<HamiltonianPtr  ,IS_CONST>::type ha_;
  typename tmptools::ConditionalAddConst<LiouvilleanPtr  ,IS_CONST>::type li_; 
  typename tmptools::ConditionalAddConst<AveragedPtr     ,IS_CONST>::type av_;
  
};



template<int RANK>
std::ostream& display(boost::shared_ptr<const Averaged<RANK> > av,
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


template<int RANK>
const LiouvilleanAveragedCommon::DArray1D average(typename LiouvilleanAveragedCommonRanked<RANK>::Ptr ptr, double t, const quantumdata::LazyDensityOperator<RANK>& matrix)
{
  return ptr ? ptr->average(t,matrix) : LiouvilleanAveragedCommon::DArray1D();
}


} // structure



#endif // STRUCTURE_STRUCTURE_H_INCLUDED
