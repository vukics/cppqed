// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines stream_densityoperator::_}
#ifndef CPPQEDCORE_QUANTUMTRAJECTORY_DENSITY_OPERATOR_STREAMER_H_INCLUDED
#define CPPQEDCORE_QUANTUMTRAJECTORY_DENSITY_OPERATOR_STREAMER_H_INCLUDED

#include "DensityOperator.h"

#include "Averaged.h"
#include "NegPT.tcc"
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
template<int RANK, typename V=tmptools::V_Empty>
class DensityOperatorStreamer
{
public:
  typedef quantumdata::DensityOperator<RANK> DensityOperator;

  typedef          structure::Averaged<RANK>      Averaged   ;
  typedef typename structure::Averaged<RANK>::Ptr AveragedPtr;

  DensityOperatorStreamer(AveragedPtr av, bool negativity) : av_(av), negativity_(negativity) {}

  StreamReturnType operator()(double t, const DensityOperator& rho, std::ostream& os, int precision) const 
  {
    auto res{structure::stream(av_,t,rho,os,precision)};
    auto & averages{std::get<1>(res)};
    if (negativity_) {
      auto n{negPT(rho,V())};
      os<<'\t'<<FormDouble(precision)(n);
#ifndef   DO_NOT_USE_FLENS
      averages.resizeAndPreserve(averages.size()+1); averages(averages.ubound(0))=n;
#endif // DO_NOT_USE_FLENS
    }
    return {os,averages};
  }

  std::ostream& streamKey(std::ostream& os, size_t& i) const
  {
    if (av_) av_->streamKey(os,i); 
    if (negativity_) os<<"Trajectory\n"<<i<<". negativity"<<std::endl;
    return os;
  }

private:
  const AveragedPtr av_ ;

  const bool negativity_;

};


} // quantumtrajectory

#endif // CPPQEDCORE_QUANTUMTRAJECTORY_DENSITY_OPERATOR_STREAMER_H_INCLUDED