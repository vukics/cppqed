// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// #include "Evolution.h"
#include "Mode.h"

#include "Master.h"
#include "QuantumJumpMonteCarlo.h"

#include <iostream>
// #include "BichromaticMode.h"
// #include "AveragingUtils.h"


const size_t cutoff=5;

using cppqedutils::json;
using namespace mode;

int main(int argc, char* argv[])
{
  // auto aOp{mode::aOp(cutoff)}, aDagOp{mode::aDagOp(cutoff)}, nOp{mode::nOp(cutoff)}, aSqrOp{aOp|aOp}, aToTheThirdOp{aOp|aSqrOp};
  //
  // std::cerr<<json(aOp)<<std::endl<<json(aDagOp)<<std::endl<<json(nOp)<<std::endl<<json(aOp|aOp|aOp)<<std::endl<<std::endl;
  //
  // auto driveOp{aOp+aDagOp};
  //
  // std::cerr<<json(driveOp)<<std::endl<<std::endl;
  //
  // {
  //   quantumdata::StateVector<1> psi({cutoff},quantumdata::zeroInit<1>), dpsidt({cutoff},quantumdata::zeroInit<1>);
  //   psi(3)={1.,1.};
  //   aOp(0.,psi,dpsidt);
  //   std::cerr<<json(dpsidt)<<std::endl<<std::endl;
  // }
  //
  // {
  //   quantumdata::StateVector<1> psi({cutoff},quantumdata::zeroInit<1>), dpsidt({cutoff},quantumdata::zeroInit<1>);
  //   psi(3)={1.,1.};
  //   aDagOp(0.,psi,dpsidt);
  //   std::cerr<<json(dpsidt)<<std::endl<<std::endl;
  // }
  //
  // {
  //   quantumdata::StateVector<3> psi({3,5,4},quantumdata::zeroInit<3>), dpsidt({3,5,4},quantumdata::zeroInit<3>);
  //
  //   psi(2,3,3)=1.;
  //
  //   for (cppqedutils::Extents<3> idx{}; idx[0]!=psi.extents[0]; cppqedutils::incrementMultiIndex(idx,psi.extents))
  //     std::cerr<<json(idx)<<" "<<psi(idx)<<std::endl;
  //
  //   quantumoperator::MultiDiagonal<3> directProduct{mode::aOp(3)*driveOp*(mode::aOp(4)|mode::aOp(4))};
  //   std::cerr<<json(directProduct)<<std::endl<<std::endl;
  //
  //   directProduct(0.,psi,dpsidt);
  //
  //   for (cppqedutils::Extents<3> idx{}; idx[0]!=psi.extents[0]; cppqedutils::incrementMultiIndex(idx,psi.extents))
  //     std::cerr<<json(idx)<<" "<<dpsidt(idx)<<std::endl;
  //
  //
  //   for (auto&& [index,offsets,diagonal] : quantumoperator::multidiagonal::range(directProduct))
  //     std::cerr/*<<json(index)*/<<json(offsets)<<json(diagonal)<<std::endl;
  // }


  auto op{optionParser()};

  cppqedutils::trajectory::Pars<quantumtrajectory::master::Pars> pt(op);

  parse(op,argc, argv);

  dcomp z{.1,1.}, eta{-.1,.2};
  double omegaKerr=2, nTh=1.;
  size_t cutoff=10;

  std::array freqs{
    structure::SystemFrequencyDescriptor{"z",z,1.},
    structure::SystemFrequencyDescriptor{"η",eta,sqrt(cutoff)},
    structure::SystemFrequencyDescriptor{"omegaKerr",omegaKerr,sqr(cutoff)}};

  cppqedutils::ODE_EngineBoost<quantumdata::StateVector<1>::StorageType> oe(quantumtrajectory::initialTimeStep(freqs),1e-12,1e-30);

  quantumtrajectory::QuantumJumpMonteCarlo
    q{QuantumSystemDynamics{
        freqs,
        std::array<structure::Lindblad<1>,2>{photonLoss(1.,nTh),photonGain(.1,nTh)},
        mode::hamiltonian(cutoff,z,omegaKerr,eta),
        expectationValues },
      quantumdata::StateVector<1>{{cutoff}},oe,randomutils::EngineWithParameters<pcg64>{1001,1},0.01};

  cppqedutils::run(q,pt,cppqedutils::trajectory::observerNoOp);


  //  for ( auto&& [idx0,idx1] : std::views::zip( cppqedutils::MultiIndexRange<5>{{5,2,4,3,2}},cppqedutils::MultiIndexRange<5>{{10,2,2,3,2}} ) )
//    std::cout<<json(idx0)<<json(idx1)<<std::endl;

  // {
  //   using namespace cppqedutils;
  //   Extents<5> ext{5,2,4,3,2};
  //   for (Extents<5> idx{}; idx[0]!=ext[0]; incrementMultiIndex(idx, ext) ) std::cout<<json(idx)<<std::endl;
  // }

//  mode::MultiDiagonal{};


/*
  // ****** Parameters of the Problem

  ParameterTable p;

  evolution::Pars<> pe(p); // Driver Parameters
  ParsBichromatic pplm(p); 

  bool
    &alternative=p.add("alternative","Alternative mode",false),
    &doStream=p.add("doStream","Stream diagonal elements of density operator",false),
    &doOffDiag=p.add("doOffDiag","Stream offdiagonal elements of density operator",false);

  // Parameter finalization
  QM_Picture& qmp=updateWithPicture(p,argc,argv);
  
  // ****** ****** ****** ****** ****** ******

  Collecting::Collection collection; collection.push_back(new AveragedQuadratures());
  if (doStream) collection.push_back(new ReducedDensityOperator<1>("",pplm.cutoff,doOffDiag));

  evolve(mode::init(pplm),(alternative ? Ptr(new DrivenDissipativeModeIP_NoExact(pplm)) : make<Collecting>(pplm,qmp,collection)),pe);
*/
}

