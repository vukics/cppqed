// -*- C++ -*-
#ifndef _TRIDIAGONAL_H
#define _TRIDIAGONAL_H

#include "TridiagonalFwd.h"

#include "DimensionsBookkeeper.h"
#include "LazyDensityOperator.h"
#include "Types.h"

#include "BlitzTinyOfArrays.h"
#include "Operators.h"

#include "Exception.h"





// Stencils???


struct TridiagonalConsistencyErrorException  : public cpputils::Exception {};
struct TridiagonalStructureMismatchException : public cpputils::Exception {};

/*
struct TridiagonalFrequenciesDiscrepancy : public TridiagonalError {};

#include "CMatrix.h"

*/


namespace quantumoperator {


template<int RANK>
inline
void
apply(const typename quantumdata::Types<RANK>::StateVectorLow& psi, typename quantumdata::Types<RANK>::StateVectorLow& dpsidt,
      const Tridiagonal<RANK>&);


const Tridiagonal<1> zero    (size_t);
const Tridiagonal<1> identity(size_t);


//////////////
//
// Tridiagonal
//
//////////////


template<int RANK> 
class Tridiagonal 
  : public DimensionsBookkeeper<RANK>, 
    private linalg::VectorSpace<Tridiagonal<RANK> >
{
public:  
  static const int LENGTH=tmptools::Power<3,RANK>::value;

  typedef blitzplusplus::TinyOfArrays<dcomp,RANK,LENGTH> Diagonals;

  typedef typename Diagonals::T_numtype Diagonal;

  typedef DimensionsBookkeeper<RANK> Base;

  typedef typename Base::Dimensions Dimensions;

  typedef typename quantumdata::Types<RANK>::StateVectorLow StateVectorLow;

  static const mpl::int_<1> _1_  ;//=mpl::int_<1>();
  static const Diagonal     empty;//=Diagonal    ();

  explicit Tridiagonal(const Diagonal& zero=empty, size_t k=0, const Diagonal& minus=empty, const Diagonal& plus=empty, mpl::int_<RANK> =_1_);

  Tridiagonal(const Tridiagonal& tridiag) 
    : Base(tridiag), diagonals_(blitzplusplus::TOA_DeepCopy(),tridiag.diagonals_), differences_(tridiag.differences_), tCurrent_(0) {}

  template<int RANK2>
  Tridiagonal(const Tridiagonal<RANK2>&, const Tridiagonal<RANK-RANK2>&); // Direct product


  template<int> // dummy template cf multiplication with Sigma
  void
  apply(const StateVectorLow& psi, StateVectorLow& dpsidt) const;


  const dcomp average(const quantumdata::LazyDensityOperator<RANK>&, 
		      const typename quantumdata::LazyDensityOperator<RANK>::Idx&, mpl::int_<RANK> =_1_);


  const Diagonals & get           () const {return diagonals_  ;}
  const Dimensions& getDifferences() const {return differences_;}


  // void           hermitianConjugateSelf();
  const Tridiagonal hermitianConjugate    () const;
  const Tridiagonal dagger                () const {return hermitianConjugate();}


  Tridiagonal& operator+=(const Tridiagonal&);

  const Tridiagonal operator-() const 
  {return Tridiagonal(*this,blitzplusplus::negate(diagonals_),differences_);}

  const Tridiagonal operator+() const {return *this;}
  Tridiagonal& operator-=(const Tridiagonal& tridiag) {(*this)+=-tridiag; return *this;}

  Tridiagonal& operator*=(const dcomp& dc);
  Tridiagonal& operator/=(const dcomp& dc) {(*this)*=1./dc; return *this;}
  Tridiagonal& operator*=(double d) {(*this)*=dcomp(d,0); return *this;}
  Tridiagonal& operator/=(double d) {(*this)*=1./dcomp(d,0); return *this;} 

  Tridiagonal& propagate(const Frequencies<RANK>&, double t);

private:
  Tridiagonal(const Base& base, const Diagonals& diagonals, const Dimensions& differences)
    : Base(base), diagonals_(blitzplusplus::TOA_DeepCopy(),diagonals), differences_(differences), tCurrent_(0) {}

  Diagonals diagonals_;

  Dimensions differences_;

  double tCurrent_;

};



template<int RANK>
inline
void
apply(const typename quantumdata::Types<RANK>::StateVectorLow& psi, typename quantumdata::Types<RANK>::StateVectorLow& dpsidt,
      const Tridiagonal<RANK>& tridiag)
{
  tridiag.apply<RANK>(psi,dpsidt);
}


template<int RANK1, int RANK2>
inline
const Tridiagonal<RANK1+RANK2>
operator*(const Tridiagonal<RANK1>& t1, const Tridiagonal<RANK2>& t2)
{
  return Tridiagonal<RANK1+RANK2>(t1,t2);
}


template<int RANK>
inline 
const Tridiagonal<RANK>
tridiagMinusHC     (const Tridiagonal<RANK>& tridiag) {return tridiag-tridiag.hermitianConjugate();}


template<int RANK>
inline 
const Tridiagonal<RANK>
tridiagPlusHC      (const Tridiagonal<RANK>& tridiag) {return tridiag+tridiag.hermitianConjugate();}


template<int RANK>
inline 
const Tridiagonal<RANK>
tridiagPlusHC_overI(const Tridiagonal<RANK>& tridiag) {return tridiagPlusHC(tridiag)/DCOMP_I;}


template<int RANK>
std::ostream& operator<<(std::ostream&, const Tridiagonal<RANK>&);


namespace details {

template<int RANK1, int RANK2, bool>
const typename Tridiagonal<RANK1+RANK2>::Diagonals
directDiagonals(const typename Tridiagonal<RANK1>::Diagonals&, const typename Tridiagonal<RANK2>::Diagonals&);

} // details


} // quantumoperator


#include "impl/Tridiagonal.tcc"


#endif // _TRIDIAGONAL_H
