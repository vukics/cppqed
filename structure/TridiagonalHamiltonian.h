// -*- C++ -*-
#ifndef _TRIDIAGONAL_HAMILTONIAN_H
#define _TRIDIAGONAL_HAMILTONIAN_H

#include "TridiagonalHamiltonianFwd.h"

#include "Hamiltonian.h"

#include "Frequencies.h"

#include<boost/tuple/tuple.hpp>

#include<list>

namespace structure {


namespace details {


template<int RANK>
class TridiagonalStore
{
protected:
  typedef quantumoperator::Tridiagonal<RANK> Tridiagonal ;
  typedef std::list<Tridiagonal>             Tridiagonals;

  TridiagonalStore(const Tridiagonals& hOverIs) : hOverIs_(hOverIs) {}

  const Tridiagonals& getH_OverIs() const {return hOverIs_;}

  mutable Tridiagonals hOverIs_;

};


template<int RANK>
class TDH_True // TridiagonalHamiltonian time-dependent base
  : public Hamiltonian<RANK,ONE_TIME>, private TridiagonalStore<RANK>
{
public:
  typedef TridiagonalStore<RANK> Base;

  typedef quantumoperator::Frequencies<RANK> Frequencies;
  typedef typename   Base::Tridiagonal       Tridiagonal;

  typedef typename Base::Tridiagonals Tridiagonals;
  typedef std::list<Frequencies>      Frequenciess;

  typedef typename Hamiltonian<RANK>::StateVectorLow StateVectorLow;

  using Base::getH_OverIs;

  TDH_True(const Tridiagonals& hOverIs, const Frequenciess& freqss) : Base(hOverIs), freqss_(freqss) {}

  virtual void addContribution(double, const StateVectorLow&, StateVectorLow&) const;

  const Frequenciess& getFreqss() const {return freqss_;}

private:
  const Frequenciess freqss_;

  using Base::hOverIs_;

};


template<int RANK>
class TDH_False // TridiagonalHamiltonian time-independent base
  : public Hamiltonian<RANK,NO_TIME>, private TridiagonalStore<RANK>
{
public:
  typedef TridiagonalStore<RANK> Base;

  typedef typename Base::Tridiagonal  Tridiagonal ;
  typedef typename Base::Tridiagonals Tridiagonals;

  typedef typename Hamiltonian<RANK>::StateVectorLow StateVectorLow;

  using Base::getH_OverIs;

  TDH_False(const Tridiagonals& hOverIs) : Base(hOverIs) {}

  virtual void addContribution(const StateVectorLow&, StateVectorLow&) const;

private:
  using Base::hOverIs_;

};


} // details


template<int RANK, bool IS_TD> // TD stands for time-dependent: the class is composed at compile-time 
class TridiagonalHamiltonian : public mpl::if_c<IS_TD,details::TDH_True<RANK>,details::TDH_False<RANK> >::type
{
public: 
  typedef typename mpl::if_c<IS_TD,details::TDH_True<RANK>,details::TDH_False<RANK> >::type Base;

  typedef typename details::TDH_True<RANK>::Frequencies  Frequencies ;
  typedef typename details::TDH_True<RANK>::Frequenciess Frequenciess;

  typedef typename Base::Tridiagonal  Tridiagonal ;
  typedef typename Base::Tridiagonals Tridiagonals;

  typedef boost::tuple<Tridiagonal,Frequencies> TridiagonalIP ;
  typedef std::list<TridiagonalIP>              TridiagonalIPs;

  // From the two sets of constructors, only one can be used, depending on the value of IS_TD
  // IS_TD=true
  TridiagonalHamiltonian(const Tridiagonal & hOverI , const Frequencies & freqs ) 
    : Base(Tridiagonals(1,hOverI),Frequenciess(1,freqs)) {}
  TridiagonalHamiltonian(const Tridiagonals& hOverIs, const Frequenciess& freqss) : Base(hOverIs,freqss) {}
  TridiagonalHamiltonian(const TridiagonalIPs&);

  // IS_TD=false
  TridiagonalHamiltonian(const Tridiagonal & hOverI                             ) 
    : Base(Tridiagonals(1,hOverI)) {}
  TridiagonalHamiltonian(const Tridiagonals& hOverIs                            ) : Base(hOverIs       ) {}

};


} // structure

#include "impl/TridiagonalHamiltonian.tcc"

#endif // _TRIDIAGONAL_HAMILTONIAN_H
