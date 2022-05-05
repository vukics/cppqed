// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines stream_densityoperator::_}
#ifndef CPPQEDCORE_QUANTUMTRAJECTORY_DENSITY_OPERATOR_STREAMER_H_INCLUDED
#define CPPQEDCORE_QUANTUMTRAJECTORY_DENSITY_OPERATOR_STREAMER_H_INCLUDED

#include "DensityOperator.h"

#include "Averaged.h"
#include "EntanglementMeasures.h"
#include "QuantumTrajectory.h"
#include "Structure.h"

#include "FormDouble.h"


namespace quantumtrajectory {


/// Wraps common functionality of Master & EnsembleMCWF concerning stream of quantum averages on the basis of density operators
/**
 * This comprises
 * - keeping a structure::Averaged instant and calling structure::Averaged::stream
 * - performing \link quantumdata::negPT negativity calculation\endlink if needed
 * - extending key with negativity when needed
 * 
 * \tparamRANK
 * \tparam V has the same function as the template parameter `V` in quantumdata::negPT
 * 
 */
template<int RANK, typename V>
class DensityOperatorStreamer
{
public:
  typedef quantumdata::DensityOperator<RANK> DensityOperator;

  DensityOperatorStreamer(::structure::AveragedPtr<RANK> av, EntanglementMeasuresSwitch ems) : av_(av), ems_(ems) {}

  StreamReturnType operator()(double t, const DensityOperator& rho, std::ostream& os, int precision) const 
  {
    auto res{structure::stream(av_,t,rho,os,precision)};
    auto & averages{std::get<1>(res)};
    if constexpr ( !isV_empty ) {
      if (ems_[0]) {
        auto n{negPT(rho,V{})};
        os<<'\t'<<FormDouble(precision)(n);
        averages.resizeAndPreserve(averages.size()+1); averages(averages.ubound(0))=n;
      }
      if (ems_[1]) {
        auto mi{mutualInformation(rho,V{})};
        os<<'\t'<<FormDouble(precision)(mi);
        averages.resizeAndPreserve(averages.size()+1); averages(averages.ubound(0))=mi;
      }
      if (ems_[2]) {
        auto p{purityOfPartialTrace(rho,V{})};
        os<<'\t'<<FormDouble(precision)(p);
        averages.resizeAndPreserve(averages.size()+1); averages(averages.ubound(0))=p;
      }

    }
    return {os,averages};
  }

  std::ostream& streamKey(std::ostream& os, size_t& i) const
  {
    if (av_) av_->streamKey(os,i); 
    if constexpr ( !isV_empty ) {
      if (ems_.any()) os<<"Trajectory\n";
      if (ems_[0]) os<<i++<<". negativity\n";
      if (ems_[1]) os<<i++<<". mutual information\n";
      if (ems_[2]) os<<i++<<". purity of partial trace\n";
    }
    return os;
  }

private:
  const ::structure::AveragedPtr<RANK> av_ ;

  const EntanglementMeasuresSwitch ems_;

  static constexpr bool isV_empty=boost::mpl::empty<V>::value;
  
};


} // quantumtrajectory

#endif // CPPQEDCORE_QUANTUMTRAJECTORY_DENSITY_OPERATOR_STREAMER_H_INCLUDED
