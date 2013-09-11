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
 * ![Inheritance diagram](../diagrams/structure.png)
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
const typename Hamiltonian<RANK>::Ptr 
qsh(DynamicsBase::Ptr base)
{return dynamic_pointer_cast<const Hamiltonian<RANK> >(base);}

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


/** \page structurebundleguide Guide on the use of the structure-bundle
 * 
 * \note The files referred to in this Section are found in directory `examples` in the distribution. The code examples are expected to compile, cf. `examples/Jamfile`.
 * 
 * ## Implementation of a class representing a harmonic-oscillator mode
 * 
 * We demonstrate how to implement an element representing a pumped lossy mode in a truncated Fock space. In the frame rotating with the pump frequency, it is described by the Hamiltonian:
 * \f[H=-\delta a^\dagger a+\lp\eta a^\dagger+\hermConj\rp,\f]
 * where \f$\delta\f$ is the detuning between the oscillator and the pump frequencies. The Liouvillean:
 * \f[\Liou\rho=\kappa\lp(n+1)\lp2a\rho a^\dagger-\comm{a^\dagger a}{\rho}_+\rp+n\lp2a^\dagger\rho a-\comm{a\,a^\dagger}{\rho}_+\rp\rp\f]
 * The frequency-like parameters are \f$\delta\f$, \f$\kappa\f$, and \f$\eta\f$, representing the mode frequency, loss rate, and pumping Rabi frequency, respectively.
 * A dimensionless parameter is \f$n\f$ (the average number of thermal photons in the heat bath) and the cutoff of the Fock space.
 * 
 * Using the notation of Sec. \ref mcwftrajectory "Description of the MCWF method" \f[J_0=\sqrt{2\kappa(n+1)}a,\quad J_1=\sqrt{2\kappa n}a^\dagger.\f]
 * The non-Hermitian Hamiltonian for the Monte Carlo wave-function method reads:
 * \f[\HnH=\lp-\delta-i\kappa(2n+1)\rp a^\dagger a+\lp\eta a^\dagger+\hermConj\rp\equiv-iz\,a^\dagger a+\lp\eta a^\dagger+\hermConj\rp,\f]
 * where we have introduced the complex frequency \f[z\equiv\kappa(2n+1)-i\delta.\f]
 * 
 * The element has to be represented by a class which inherits publicly from the necessary classes in the structure namespace.
 * In this simple case, it is basically two helper functions returning quantumoperator::Tridiagonal instances, a constructor, and two virtual functions inherited from
 * structure::ElementAveraged that have to be written.
 * 
 * \see The description of quantumoperator::Tridiagonal
 * 
 * Consider the file `ExampleMode.h`:
 * 
 * \snippet ExampleMode.h basic example mode
 * 
 * This will suffice here. Let us look at the implementations in `ExampleMode.cc`:\dontinclude ExampleMode.cc
 * \until sqrt(cutoff))),
 * We construct the Free base with the dimension of the system and the name-value-multiplier tuples for the frequency-like parameters of the system, which in this case are all complex
 * (cf. the explanation @ DynamicsBase). #FREQS is a macro defined in DynamicsBase.h, which facilitates the creation of lists of such tuples.
 * \until cutoff))),
 * We construct the time-independent TridiagonalHamiltonian base. This is greatly facilitated by the algebra and helpers of the quantumoperator::Tridiagonal class.
 * \warning When implementing the Hamiltonian, not \f$H\f$ itself but \f$\frac Hi\f$ has to supplied!
 * 
 * \until ("photon absorption")),
 * We construct the ElementLiouvillean base whose second template argument denotes the number of different quantum jumps, which is 2 in this case.
 * The constructor takes the strategies for calculating the impact of a jump on a free::StateVectorLow, and for calculating the rate from a free::LazyDensityOperator.
 * These strategy functions are produced from the free-standing helpers in Lines 10-14 above through argument binding.
 * \until ("imag(\")"))
 * We construct the ElementAveraged base, with parameters necessary to produce a simple key for quantum averages communicated towards the user.
 * Here we calculate only three such averages, the expectation value of the number operator, and the real and imaginary parts of that of the ladder operator.
 * \until }
 * With the DynamicsBase::getParsStream function we obtain a stream whereon we can write more information about the object that gets communicated towards the user
 * in that part of the output which summarizes the parameters of the actual run.
 *
 * Next, the inherited function Averaged::average_v is implemented:
 * \skip average_v
 * \until aJumpRate
 * the expectation value of the photon number is calculated (where we can reuse our function `aJumpRate`, with unit loss rate).
 * \until imag(offdiag);
 * the expectation value of the ladder operator is calculated (real & imaginary parts)
 * \until }
 * \until }
 * 
 * The implementation of the helpers is also quite straightforward. It may come to a separate file `ExampleModeImpl.cc`:\dontinclude ExampleModeImpl.cc
 * \until blitz::tensor::i+1.
 * \warning In situations like this, it is extremely important to write `blitz::tensor::i+1.`, otherwise the `blitz::sqrt` function will operate within the integers,
 * which is perhaps a shortcoming of the Blitz++ library in this respect.
 * 
 * \until const Tridiagonal nop
 * \until }
 * 
 * ### Exploiting interaction picture
 * 
 * In many situations, it pays to transfer to interaction picture defined by the first term of the Hamiltonian \f$\HnH\f$ above. The ladder operator in interaction picture is
 * \f[a\Int(t)=a e^{-zt},\f] so that the Hamiltonian reads \f[H\Int(t)=\lp\eta a^\dagger e^{zt}+\eta^*ae^{-zt}\rp.\f]
 * 
 * \see [These notes](http://optics.szfki.kfki.hu/~vukics/Pictures.pdf) on how to treat interaction pictures defined by non-unitary transition operators in a consistent way.
 * 
 * In this case, the class representing the element has to be derived from Exact as well, which represents the transformation between the two pictures.
 * In addition, instead of TridiagonalHamiltonian`<1,false>`, we need to derive from TridiagonalHamiltonian`<1,true>` because the Hamiltonian is now time-dependent.
 * \note In general usage, the jump and the averages are calculated in the normal picture also in this case (cf. explanation of classes Hamiltonian, Exact, Liouvillean, and Averaged,
 * and furthermore quantumtrajectory::MCWF_Trajectory). This allows for reusing the same code in both pictures. (Incidentally, here the Liouvillean remains unchanged anyway.)
 * 
 * \dontinclude ExampleMode.h
 * \skip FreeExact
 * \until PumpedLossyModeIP&);
 * In the implementation, the only difference from the previous case will be the constructor, because the Hamiltonian now also requires furnishing with frequencies
 * (cf. quantumoperator::furnishWithFreqs), and the implementation of the virtual function FreeExact::updateU.
 * 
 * When “furnished with frequencies”, a quantumoperator::Tridiagonal object will internally take care about the time-dependent phases appearing in \f$a\Int(t)\f$ and \f$H\Int(t)\f$:
 * \dontinclude ExampleMode.cc
 * \skip mainDiagonal
 * \until blitz::tensor::i));
 * \until }
 * FreeExact assumes that the operator transforming between the two pictures is diagonal, and the factors to update are simply its diagonal elements.
 * If this is not the case, Exact has to be used instead.
 * 
 * The construction of the object representing the ladder operator furnished with frequencies is also straightforward:\dontinclude ExampleModeImpl.cc
 * \skip aop(const PumpedLossyModeIP& mode)
 * \until }
 * \note Since a lot of the code from the previous case can be reused here, one will usually adopt an inheritence- or class-composition-based solution to implement classes like
 * `PumpedLossyMode` and `PumpedLossyModeIP` (cf. the actual implementation of a harmonic-oscillator mode in the framework in `elements/frees/Mode.h`).
 * 
 * ## Implementing an X-X interaction
 * 
 * Let us consider the interaction described by the Hamiltonian \f[H_\text{X-X}=g(a+a^\dagger)(b+b^\dagger).\f]
 * 
 * The class implementing this interaction has to be derived from Interaction`<2>` because it is a binary interaction, and TridiagonalHamiltonian`<2,...>`
 * (note that quantumoperator::Tridiagonal is capable to represent direct products of tridiagonal matrices).
 * 
 * The only thing requiring some care is that once we transform some elements into interaction picture, the whole Hamiltonian is transformed, that is,
 * \f$a\f$ or \f$b\f$ or both may be in interaction picture. Here, for the sake of simplicity, we assume that both constituents are of the type `PumpedLossyModeIP`.
 * Hence, the Hamiltonian is in fact \f[H_{\text{X-X;I}}(t)=g(ae^{-z_at}+a^\dagger e^{z_at})(be^{-z_bt}+b^\dagger e^{z_bt}).\f]
 * 
 * Consider `ExampleInteraction.h`:\dontinclude ExampleInteraction.h
 * \skip #include
 * \until };
 * 
 * `ExampleInteraction.cc` then reads \include ExampleInteraction.cc
 * As we see, the Hamiltonian can be written in a rather straightforward way, and it internally takes care about the time-dependent phases appearing in \f$H_{\text{X-X;I}}(t)\f$,
 * which result from the use of interaction picture.
 *
 */


} // structure



#endif // STRUCTURE_STRUCTURE_H_INCLUDED
