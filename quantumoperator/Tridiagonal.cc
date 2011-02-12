// #include<boost/mpl/int.hpp>

#include "Tridiagonal.h"


// namespace mpl=boost::mpl;


namespace quantumoperator {

template<>
const mpl::int_<1> Tridiagonal<1>::_1_ =mpl::int_<1>();

template<>
const Tridiagonal<1>::Diagonal Tridiagonal<1>::empty=Tridiagonal<1>::Diagonal();


/////////
//
// RANK=1
//
/////////

template<> template<>
void Tridiagonal<1>::apply<1>(const StateVectorLow& psi, StateVectorLow& dpsidt) const
{
  using blitz::Range;

  const Diagonals & a=diagonals_  ;
  const Dimensions& K=differences_;

  Range 
    I0l(0,psi.ubound(0)-int(K(0))),
    I0h(I0l+int(K(0)));


  if (a(0).size()) dpsidt     +=a(0)     *psi     ;
  if (a(1).size()) dpsidt(I0h)+=a(1)(I0l)*psi(I0l);
  if (a(2).size()) dpsidt(I0l)+=a(2)(I0l)*psi(I0h);

  // NEEDS_WORK consider whether    
  // dpsidt(I0) += a(0)(I0     )*psi(I0     ) + a(1)(I0-K(0))*psi(I0-K(0)) + a(2)(I0     )*psi(I0+K(0));
  // can be made work with some conditional expression (a dummy array expression?)
  
}


/////////
//
// RANK=2
//
/////////

template<> template<>
void Tridiagonal<2>::apply<2>(const StateVectorLow& psi, StateVectorLow& dpsidt) const
{
  using blitz::Range;
  
  const Diagonals & a=diagonals_  ;
  const Dimensions& K=differences_;

  Range 
    I0(0,psi.ubound(0)),
    I0l(0,psi.ubound(0)-int(K(0))),
    I0h(I0l+int(K(0))),
    I1(0,psi.ubound(1)),
    I1l(0,psi.ubound(1)-int(K(1))),
    I1h(I1l+int(K(1)));


  if (a(0).size()) dpsidt         +=a(0)         *psi         ;
  if (a(1).size()) dpsidt(I0 ,I1h)+=a(1)(I0 ,I1l)*psi(I0 ,I1l);
  if (a(2).size()) dpsidt(I0 ,I1l)+=a(2)(I0 ,I1l)*psi(I0 ,I1h);
  if (a(3).size()) dpsidt(I0h,I1 )+=a(3)(I0l,I1 )*psi(I0l,I1 );
  if (a(4).size()) dpsidt(I0h,I1h)+=a(4)(I0l,I1l)*psi(I0l,I1l);
  if (a(5).size()) dpsidt(I0h,I1l)+=a(5)(I0l,I1l)*psi(I0l,I1h);
  if (a(6).size()) dpsidt(I0l,I1 )+=a(6)(I0l,I1 )*psi(I0h,I1 );
  if (a(7).size()) dpsidt(I0l,I1h)+=a(7)(I0l,I1l)*psi(I0h,I1l);
  if (a(8).size()) dpsidt(I0l,I1l)+=a(8)(I0l,I1l)*psi(I0h,I1h);

}



/////////
//
// RANK=3
//
/////////

template<> template<>
void Tridiagonal<3>::apply<3>(const StateVectorLow& psi, StateVectorLow& dpsidt) const
{
  using blitz::Range;
  
  const Diagonals & a=diagonals_  ;
  const Dimensions& K=differences_;

  Range 
    I0(0,psi.ubound(0)),
    I0l(0,psi.ubound(0)-int(K(0))),
    I0h(I0l+int(K(0))),
    I1(0,psi.ubound(1)),
    I1l(0,psi.ubound(1)-int(K(1))),
    I1h(I1l+int(K(1))),
    I2(0,psi.ubound(2)),
    I2l(0,psi.ubound(2)-int(K(2))),
    I2h(I2l+int(K(2)));


  if (a(0 ).size()) dpsidt             +=a(0 )             *psi             ;
  if (a(1 ).size()) dpsidt(I0 ,I1 ,I2h)+=a(1 )(I0 ,I1 ,I2l)*psi(I0 ,I1 ,I2l);
  if (a(2 ).size()) dpsidt(I0 ,I1 ,I2l)+=a(2 )(I0 ,I1 ,I2l)*psi(I0 ,I1 ,I2h);
  if (a(3 ).size()) dpsidt(I0 ,I1h,I2 )+=a(3 )(I0 ,I1l,I2 )*psi(I0 ,I1l,I2 );
  if (a(4 ).size()) dpsidt(I0 ,I1h,I2h)+=a(4 )(I0 ,I1l,I2l)*psi(I0 ,I1l,I2l);
  if (a(5 ).size()) dpsidt(I0 ,I1h,I2l)+=a(5 )(I0 ,I1l,I2l)*psi(I0 ,I1l,I2h);
  if (a(6 ).size()) dpsidt(I0 ,I1l,I2 )+=a(6 )(I0 ,I1l,I2 )*psi(I0 ,I1h,I2 );
  if (a(7 ).size()) dpsidt(I0 ,I1l,I2h)+=a(7 )(I0 ,I1l,I2l)*psi(I0 ,I1h,I2l);
  if (a(8 ).size()) dpsidt(I0 ,I1l,I2l)+=a(8 )(I0 ,I1l,I2l)*psi(I0 ,I1h,I2h);
  if (a(9 ).size()) dpsidt(I0h,I1 ,I2 )+=a(9 )(I0l,I1 ,I2 )*psi(I0l,I1 ,I2 );
  if (a(10).size()) dpsidt(I0h,I1 ,I2h)+=a(10)(I0l,I1 ,I2l)*psi(I0l,I1 ,I2l);
  if (a(11).size()) dpsidt(I0h,I1 ,I2l)+=a(11)(I0l,I1 ,I2l)*psi(I0l,I1 ,I2h);
  if (a(12).size()) dpsidt(I0h,I1h,I2 )+=a(12)(I0l,I1l,I2 )*psi(I0l,I1l,I2 );
  if (a(13).size()) dpsidt(I0h,I1h,I2h)+=a(13)(I0l,I1l,I2l)*psi(I0l,I1l,I2l);
  if (a(14).size()) dpsidt(I0h,I1h,I2l)+=a(14)(I0l,I1l,I2l)*psi(I0l,I1l,I2h);
  if (a(15).size()) dpsidt(I0h,I1l,I2 )+=a(15)(I0l,I1l,I2 )*psi(I0l,I1h,I2 );
  if (a(16).size()) dpsidt(I0h,I1l,I2h)+=a(16)(I0l,I1l,I2l)*psi(I0l,I1h,I2l);
  if (a(17).size()) dpsidt(I0h,I1l,I2l)+=a(17)(I0l,I1l,I2l)*psi(I0l,I1h,I2h);
  if (a(18).size()) dpsidt(I0l,I1 ,I2 )+=a(18)(I0l,I1 ,I2 )*psi(I0h,I1 ,I2 );
  if (a(19).size()) dpsidt(I0l,I1 ,I2h)+=a(19)(I0l,I1 ,I2l)*psi(I0h,I1 ,I2l);
  if (a(20).size()) dpsidt(I0l,I1 ,I2l)+=a(20)(I0l,I1 ,I2l)*psi(I0h,I1 ,I2h);
  if (a(21).size()) dpsidt(I0l,I1h,I2 )+=a(21)(I0l,I1l,I2 )*psi(I0h,I1l,I2 );
  if (a(22).size()) dpsidt(I0l,I1h,I2h)+=a(22)(I0l,I1l,I2l)*psi(I0h,I1l,I2l);
  if (a(23).size()) dpsidt(I0l,I1h,I2l)+=a(23)(I0l,I1l,I2l)*psi(I0h,I1l,I2h);
  if (a(24).size()) dpsidt(I0l,I1l,I2 )+=a(24)(I0l,I1l,I2 )*psi(I0h,I1h,I2 );
  if (a(25).size()) dpsidt(I0l,I1l,I2h)+=a(25)(I0l,I1l,I2l)*psi(I0h,I1h,I2l);
  if (a(26).size()) dpsidt(I0l,I1l,I2l)+=a(26)(I0l,I1l,I2l)*psi(I0h,I1h,I2h);

}



namespace {

bool Compatible(const TTD_CARRAY(1)& a, size_t diffA, const TTD_CARRAY(1)& b, size_t diffB=0)
{
  if (!a.size() || !b.size() || a.size()-diffA==b.size()-diffB) return true;
  return false;
}

} 


template<>
Tridiagonal<1>::Tridiagonal(const Diagonal& zero, size_t k, const Diagonal& minus, const Diagonal& plus, mpl::int_<1>)
  : Base(std::max(std::max(size_t(zero.size()),minus.size()+k),plus.size()+k)),
    // funnily enough, Array::size() returns an int ...
    diagonals_(blitzplusplus::TOA_DeepCopy(),zero,minus,plus),
    differences_(k),
    tCurrent_(0)
{
  // Consistency check:
  if (
      !Compatible(minus,k,zero) || !Compatible(plus,k,zero) || !Compatible(minus,k,plus,k) || 
      ((minus.size() || plus.size()) && !k)
      )
    throw TridiagonalConsistencyErrorException();

}


const Tridiagonal<1> identity(size_t dim)
{
  Tridiagonal<1>::Diagonal diagonal(dim);
  diagonal=1.;
  return Tridiagonal<1>(diagonal);
}



namespace details {


size_t binOp1(size_t otherDifference, size_t& difference)
// Declared in Tridiagonal.tcc
{
  if (!difference) return otherDifference;
  else if (difference==otherDifference) return difference;
  else throw TridiagonalStructureMismatchException();
}


} // details




} // quantumoperator


