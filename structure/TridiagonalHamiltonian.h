// -*- C++ -*-
#ifndef _TRIDIAGONAL_HAMILTONIAN_H
#define _TRIDIAGONAL_HAMILTONIAN_H

#include "TridiagonalHamiltonianFwd.h"

#include "Hamiltonian.h"

#include "Tridiagonal.h"

#include <list>

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
        Tridiagonals& getH_OverIs()       {return const_cast<Tridiagonals&>(static_cast<const TridiagonalStore*>(this)->getH_OverIs());}

  mutable Tridiagonals hOverIs_;

};


template<int RANK>
class TDH_True // TridiagonalHamiltonian time-dependent base
  : public Hamiltonian<RANK,ONE_TIME>, private TridiagonalStore<RANK>
{
public:
  typedef TridiagonalStore<RANK> Base;

  typedef typename Base::Tridiagonal  Tridiagonal;
  typedef typename Base::Tridiagonals Tridiagonals;

  typedef typename Hamiltonian<RANK>::StateVectorLow StateVectorLow;

  using Base::getH_OverIs;

  TDH_True(const Tridiagonals& hOverIs) : Base(hOverIs) {}

  virtual void addContribution(double, const StateVectorLow&, StateVectorLow&) const;

private:
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

  typedef typename Base::Tridiagonal  Tridiagonal ;
  typedef typename Base::Tridiagonals Tridiagonals;

  TridiagonalHamiltonian(const Tridiagonal & hOverI ) : Base(Tridiagonals(1,hOverI)) {}
  TridiagonalHamiltonian(const Tridiagonals& hOverIs) : Base(hOverIs) {}

};


} // structure


#endif // _TRIDIAGONAL_HAMILTONIAN_H
