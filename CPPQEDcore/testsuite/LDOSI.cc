// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "LazyDensityOperatorSliceIterator.h"
#include "Mode_.h"

static const int R=2;

typedef quantumdata::StateVector        <R> StateVector        ;
typedef quantumdata::DensityOperator    <R> DensityOperator    ;
typedef quantumdata::LazyDensityOperator<R> LazyDensityOperator;

typedef DensityOperator::Idx Idx;

template<typename V>
struct LDOSI_MF : mpl::identity<quantumdata::ldo::DiagonalIterator<R,V> > {};

template<typename V>
struct LDO_Res_MF : mpl::identity<quantumdata::LazyDensityOperator<mpl::size<V>::value> > {};

const tmptools::Vector<0> vec0=tmptools::Vector<0>();
const tmptools::Vector<1> vec1=tmptools::Vector<1>();


using namespace std;

ostream& operator<<(ostream& os, const quantumdata::LazyDensityOperator<1>& ldo)
{
  const size_t dim=ldo.getDimension();
  for (size_t i=0; i<dim; ++i) {
    for (size_t j=0; j<dim; ++j) 
      os<<ldo(i,j)<<'\t';
    os<<endl;
  }
  return os;
}

void print(const quantumdata::LazyDensityOperator<1>& ldo)
{
  cout<<ldo<<endl;
}

int main()
{
  StateVector psi(StateVector::Dimensions(5,3));
  psi()(2,1)=dcomp(1, 2);
  psi()(2,2)=dcomp(2,-2);
  psi()(4,1)=dcomp(3, 1);
  psi()(4,2)=dcomp(1,-2);
  
  // LDOSI_MF<Vec1>::type i(psi.begin<Vec1>());

  // const LDO_Res_MF<Vec1>::type& ldo(*i);

  // cout<<*++++i;

  std::for_each(psi.begin(vec0),psi.end(vec0),bind(print,_1));
  std::for_each(psi.begin(vec1),psi.end(vec1),bind(print,_1));

  DensityOperator rho(StateVector::Dimensions(3,4));
  rho()(Idx(2,1,2,2))=dcomp(1,2);
  rho()(Idx(1,3,0,3))=dcomp(1,0);

  std::for_each(rho.begin(vec0),rho.end(vec0),bind(print,_1));
  std::for_each(rho.begin(vec1),rho.end(vec1),bind(print,_1));

}
