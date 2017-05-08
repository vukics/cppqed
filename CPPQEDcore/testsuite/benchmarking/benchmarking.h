// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Composite.h"

#include "Randomized.h"

#include <boost/progress.hpp>


using namespace std;
using namespace randomized;
using namespace tmptools;


template<int RANK, typename V, typename SYS>
void benchmark(const structure::QuantumSystem<RANK>& sys, const SYS& sys_v, V v, bool doDisplay=false, size_t nRepeat=1000)
{
  static const int RANK_V=mpl::size<V>::type::value;
  typedef structure::Hamiltonian<RANK_V> Ha_V;

  sys.displayParameters(std::cout);

  quantumdata::StateVector<RANK> psi(sys.getDimensions());

  {
    Randomized::Ptr Ran(MakerGSL()(1001));
    boost::generate(psi(),bind(&Randomized::dcompRan,Ran));
  }

  // 1

  if (doDisplay) cout<<blitzplusplus::unaryArray(psi());
  quantumdata::StateVector<RANK> psiout(psi);

  {
    const structure::Hamiltonian<RANK>* ha=dynamic_cast<const structure::Hamiltonian<RANK>*>(&sys);

    boost::progress_timer t;
    for (unsigned count=0; count<nRepeat; ++count)
      structure::Hamiltonian<RANK>::addContribution(0.,psi(),psiout(),0.,ha);
  }

  if (doDisplay) cout<<blitzplusplus::unaryArray(psiout());


  // 2

  const Ha_V* ha=dynamic_cast<const Ha_V*>(&sys_v);

  {
    // quantumdata::StateVector<RANK> psiout(psi,false);

    boost::progress_timer t;
    for (unsigned count=0; count<nRepeat; ++count)
      cpputils::for_each(blitzplusplus::basi::fullRange(psi(),v),blitzplusplus::basi::begin(psiout(),v),
                         bind(&Ha_V::addContribution,ha,0,_1,_2,0));
  }

  // 3

  {
    const blitzplusplus::SlicesData<RANK,V> sd(psi());

    boost::progress_timer t;
    for (unsigned count=0; count<nRepeat; ++count)
      cpputils::for_each(blitzplusplus::basi_fast::fullRange(sd,psi(),v),blitzplusplus::basi_fast::begin(sd,psiout(),v),
                         bind(&Ha_V::addContribution,ha,0,_1,_2,0));
  }


}
