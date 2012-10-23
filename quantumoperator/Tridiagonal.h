// -*- C++ -*-
#ifndef QUANTUMOPERATOR_TRIDIAGONAL_H_INCLUDED
#define QUANTUMOPERATOR_TRIDIAGONAL_H_INCLUDED


#ifndef   QUANTUMOPERATOR_TRIDIAGONAL_MAX_RANK
#define   QUANTUMOPERATOR_TRIDIAGONAL_MAX_RANK BLITZ_ARRAY_LARGEST_RANK
#endif // QUANTUMOPERATOR_TRIDIAGONAL_MAX_RANK

#include "TridiagonalFwd.h"

#include "DimensionsBookkeeper.h"
#include "LazyDensityOperator.h"
#include "Types.h"

#include "impl/BlitzTinyOfArrays.tcc"
#include "Operators.h"

#include "Exception.h"





// Stencils???


struct TridiagonalConsistencyErrorException  : public cpputils::Exception {};
struct TridiagonalStructureMismatchException : public cpputils::Exception {};
struct TridiagonalTimeMismatchException      : public cpputils::Exception {};

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


template<int RANK>
const Tridiagonal<RANK>
furnishWithFreqs(const Tridiagonal<RANK>& tridiag, const typename Tridiagonal<RANK>::Diagonal& mainDiagonal);


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
  static const int N_RANK=RANK;

  typedef blitzplusplus::TinyOfArrays<dcomp,RANK,LENGTH> Diagonals;

  typedef typename Diagonals::T_numtype Diagonal;

  typedef DimensionsBookkeeper<RANK> Base;

  typedef typename Base::Dimensions Dimensions;

  typedef typename quantumdata::Types<RANK>::StateVectorLow StateVectorLow;

private:
  typedef mpl::int_<RANK> IntRANK;
  typedef blitz::TinyVector<blitz::TinyVector<blitz::Range,3>,RANK> Ranges;

  static const mpl::int_<1> _1_  ;//=mpl::int_<1>();
  static const Diagonal     empty;//=Diagonal    ();

public:
  explicit Tridiagonal(const Diagonal& zero=empty, size_t k=0, const Diagonal& minus=empty, const Diagonal& plus=empty, bool toFreqs=false, IntRANK=_1_);

  Tridiagonal(const Tridiagonal& tridiag) 
    : Base(tridiag), diagonals_(blitzplusplus::TOA_DeepCopy(),tridiag.diagonals_), differences_(tridiag.differences_), tCurrent_(tridiag.tCurrent_), 
      freqs_(blitzplusplus::TOA_DeepCopy(),tridiag.freqs_) {}

  template<int RANK2>
  Tridiagonal(const Tridiagonal<RANK2>&, const Tridiagonal<RANK-RANK2>&); // Direct product


  void apply(const StateVectorLow& psi, StateVectorLow& dpsidt) const;


  const dcomp average(const quantumdata::LazyDensityOperator<RANK>&, 
		      const typename quantumdata::LazyDensityOperator<RANK>::Idx&, IntRANK=_1_);


  const Diagonals & get           () const {return   diagonals_;}
  const Dimensions& getDifferences() const {return differences_;}

  double            getTime       () const {return    tCurrent_;}

  const Diagonals & getFreqs      () const {return       freqs_;}


  // void           hermitianConjugateSelf();
  const Tridiagonal hermitianConjugate    () const;
  const Tridiagonal dagger                () const {return hermitianConjugate();}

  Tridiagonal& furnishWithFreqs(const Diagonal& mainDiagonal, IntRANK=_1_);

  Tridiagonal& operator+=(const Tridiagonal&);

  const Tridiagonal operator-() const 
  {return Tridiagonal(*this,blitzplusplus::negate(diagonals_),differences_,tCurrent_,freqs_);}

  const Tridiagonal operator+() const {return *this;}
  Tridiagonal& operator-=(const Tridiagonal& tridiag) {(*this)+=-tridiag; return *this;}

  Tridiagonal& operator*=(const dcomp& dc);
  Tridiagonal& operator/=(const dcomp& dc) {(*this)*=1./dc; return *this;}
  Tridiagonal& operator*=(double d) {(*this)*=dcomp(d,0); return *this;}
  Tridiagonal& operator/=(double d) {(*this)*=1./dcomp(d,0); return *this;} 

  Tridiagonal& propagate(double t);

private:
  Tridiagonal(const Base& base, const Diagonals& diagonals, const Dimensions& differences, double tCurrent, const Diagonals& freqs)
    : Base(base), diagonals_(blitzplusplus::TOA_DeepCopy(),diagonals), differences_(differences), tCurrent_(tCurrent), freqs_(blitzplusplus::TOA_DeepCopy(),freqs) {}

  ///////////////////////
  // Apply implementation
  ///////////////////////

  template<int START, typename V_DPSIDT, typename V_A, typename V_PSI, int REMAINING>
  void doApply(mpl::int_<REMAINING>,const Ranges&, const StateVectorLow&, StateVectorLow&) const;

  // "Specialization" for the REMAINING=0 case:
  // This should be specialized for all RANKs in Tridiagonal.tcc
  template<int START, typename V_DPSIDT, typename V_A, typename V_PSI>
  void doApply(mpl::int_<0>,const Ranges&, const StateVectorLow&, StateVectorLow&) const;

  struct FillRangesHelper
  {
    typedef const typename StateVectorLow::T_index Bound;
    
    FillRangesHelper(Ranges& ranges, const Bound& ubound, const Dimensions& k) : ranges_(ranges), ubound_(ubound), k_(k) {}

    template<typename ICW> void operator()(ICW);
  
  private:
    Ranges& ranges_;
    const Bound& ubound_;
    const Dimensions& k_;

  };

  const Ranges fillRanges(const typename StateVectorLow::T_index&) const;

  ///////
  // Data
  ///////

  Diagonals diagonals_;

  Dimensions differences_;

  double tCurrent_;

  Diagonals freqs_;

};



template<int RANK>
inline
void
apply(const typename quantumdata::Types<RANK>::StateVectorLow& psi, typename quantumdata::Types<RANK>::StateVectorLow& dpsidt,
      const Tridiagonal<RANK>& tridiag)
{
  tridiag.apply(psi,dpsidt);
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


#endif // QUANTUMOPERATOR_TRIDIAGONAL_H_INCLUDED
