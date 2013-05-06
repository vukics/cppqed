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
const typename Hamiltonian<RANK,TWO_TIME>::Ptr 
qsh(boost::shared_ptr<const QuantumSystem<RANK> > quantumSystem)
{return dynamic_pointer_cast<const Hamiltonian<RANK,TWO_TIME> >(quantumSystem);}

template<int RANK>
inline 
const typename Liouvillean<RANK,true>::Ptr 
qsl(boost::shared_ptr<const QuantumSystem<RANK> > quantumSystem)
{return dynamic_pointer_cast<const Liouvillean<RANK,true> >(quantumSystem);}

template<int RANK>
inline 
const typename Averaged<RANK,true>::Ptr 
qsa(boost::shared_ptr<const QuantumSystem<RANK> > quantumSystem)
{return dynamic_pointer_cast<const Averaged<RANK,true> >(quantumSystem);}



template<int RANK>
inline 
const typename Exact<RANK>::Ptr 
qse(DynamicsBase::Ptr base)
{return dynamic_pointer_cast<const Exact<RANK> >(base);}

template<int RANK>
inline 
const typename Hamiltonian<RANK,TWO_TIME>::Ptr 
qsh(DynamicsBase::Ptr base)
{return dynamic_pointer_cast<const Hamiltonian<RANK,TWO_TIME> >(base);}

template<int RANK>
inline 
const typename Liouvillean<RANK,true>::Ptr 
qsl(DynamicsBase::Ptr base)
{return dynamic_pointer_cast<const Liouvillean<RANK,true> >(base);}

template<int RANK>
inline 
const typename Averaged<RANK,true>::Ptr 
qsa(DynamicsBase::Ptr base)
{return dynamic_pointer_cast<const Averaged<RANK,true> >(base);}


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
  typedef typename Ex::Ptr ExactPtr;
  typedef typename Ha::Ptr HamiltonianPtr;
  typedef typename Li::Ptr LiouvilleanPtr;
  typedef typename Av::Ptr AveragedPtr;

  typedef typename Ex::StateVectorLow StateVectorLow;

  typedef typename Li::Probabilities       Probabilities      ;
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

  void actWithU(double t, StateVectorLow& psi) const {if (ex_) ex_->actWithU(t,psi);}


  // Hamiltonian
  
  void addContribution(double t, const StateVectorLow& psi, StateVectorLow& dpsidt, double tIntPic0) const {if (ha_) ha_->addContribution(t,psi,dpsidt,tIntPic0);}


  // Liouvillean
  
  void actWithJ(double t, StateVectorLow& psi, size_t jumpNo) const {if (li_) li_->actWithJ(t,psi,jumpNo);}

  
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
