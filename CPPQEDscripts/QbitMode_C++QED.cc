// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "JaynesCummings.h"
#include "Qbit.h"

#include "QuantumJumpMonteCarlo.h"


using namespace std; using namespace quantumtrajectory;


int main(int argc, char* argv[])
{
  auto op{optionParser()};
  
  cppqedutils::trajectory::Pars<quantumtrajectory::qjmc::Pars<pcg64>> pt(op);
  
  qbit::Pars pq(op,"Q");
  mode::Pars pm(op,"M");
  
  parse(op,argc, argv);

  BinarySystem bs{
    make(pq),2,make(pm),pm.cutoff,
    SystemFrequencyStore{},Liouvillian<2>{},
    structure::makeHamiltonianElement<2>("jc",jaynescummings::hamiltonian<0,1>(pm.cutoff,1.)),
    exact_propagator_ns::noOp,
    expectation_values_ns::noOp};


  quantumdata::StateVector<2> psi{{2,pm.cutoff}}; psi(0,0)=1;// psi(1)=1; 
 
  run(
    qjmc::make<cppqedutils::ODE_EngineBoost>(bs,/*quantumdata::StateVector<1>{{pm.cutoff}}*/std::move(psi),pt),
    // quantumtrajectory::qjmc::makeEnsemble<cppqedutils::ODE_EngineBoost>(std::move(mode),/*quantumdata::StateVector<1>{{pm.cutoff}}*/std::move(psi),pt),
    pt,trajectory::observerNoOp);
  
  
  // evolution::Pars<> pe(p); // Driver Parameters
  // qbit::ParsDrivenDissipativePhaseNoise pplqb(p); 
  // mode::ParsDrivenDissipative pplm (p); 
  // jaynescummings::Pars  pjc  (p); 
  // 
  // // Parameter finalization
  // QM_Picture& qmp=updateWithPicture(p,argc,argv);
  
  // ****** ****** ****** ****** ****** ******

  /*
  mode::StateVector psiMode(pplm.cutoff);
  psiMode=(mode::fock(0,pplm.cutoff,mathutils::PI/2.)+mode::fock(1,pplm.cutoff)+mode::fock(2,pplm.cutoff,-mathutils::PI/2.))/sqrt(3.);
  */
  // entangled state: 
  // psi=qbit::state0()*mode::fock(1,pplm.cutoff)+qbit::state1()*mode::fock(0,pplm.cutoff);
  // psi.renorm();

  // evolve(qbit::init(pplqb)*mode::init(pplm),
  //        binary::make(jaynescummings::make(qbit::make(pplqb,qmp),
  //                                          mode::make(pplm ,qmp/*,structure::averaged::ReducedDensityOperator("Mode",pplm.cutoff,true)*/),
  //                                          pjc)),
  //        pe);


}



/*
void f()
{
  static_assert( std::forward_iterator<SliceIterator<std::span<const size_t>,dcomp,3>> );
  static_assert( std::forward_iterator<SliceIterator<std::vector<size_t>,dcomp,3>> );
  
  static_assert( std::ranges::range<SliceRangeReferencing<dcomp,3>> );
  static_assert( std::ranges::range<SliceRangeOwning<dcomp,3>> );
  
  static_assert( std::ranges::viewable_range<SliceRangeOwning<dcomp,3>> );

  std::vector<size_t> offsets{0,10,12,20,22,30,32,40};
  
  MultiArray<dcomp,5> array{{2,4,2,3,5}}, array2{{2,4,2,3,5}};
  
  auto sr=sliceRange<retainedAxes<1,3,0>>(array), sr2=sliceRange<retainedAxes<1,3,0>>(array2);
  
  auto probe = std::views::all(sliceRange<retainedAxes<1,3,0>>(array,offsets));
  
  auto probe2 = std::views::zip(sr,sr2);
  
  auto probe3 = std::views::zip(sliceRange<retainedAxes<1,3,0>>(array),sliceRange<retainedAxes<1,3,0>>(array2));

  std::vector<size_t> offsets{0,10,12,20,22,30,32,40};
  
  MultiArray<dcomp,5> array{{2,4,2,3,5}}, array2{{2,4,2,3,5}};
  
  auto sr=sliceRange<retainedAxes<2,1,3,0,4>>(array)//, sr2=sliceRange<retainedAxes<1,3,0>>(array2)
  ;

  // Print the extent
  std::cout << "The extent of the range is: " << std::ranges::distance(std::ranges::begin(sr),std::ranges::end(sr)) << std::endl;

}
*/


//   const std::vector<size_t> offsets0;
//   
// //  static_assert( expectation_values_ns::functional< decltype(getEV(bs)), 2 > );
//     
//   quantumdata::StateVector<2> psi{{2,pm.cutoff}};// psi(0)=0; psi(1)=1;
//   
//   partialTrace<retainedAxes<0>,2>(
//     StateVectorConstView<2>(psi.constView()),
//     offsets0,
//     [&] (auto m) {return expectation_values_ns::calculate<1>(getEV(bs.qsd0),0.,m);},
//     [] (const auto& a, const auto& b) {return a;}
//   );

