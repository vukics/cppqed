#include "Tridiagonal.h"
#include "StateVector.h"

#include "Randomized.h"

#include <boost/progress.hpp>

#include <boost/preprocessor/arithmetic/sub.hpp>
#include <boost/preprocessor/comparison/less.hpp>
#include <boost/preprocessor/control/iif.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>

#define RANK 5

using namespace std;

const double totalDim=59049;// 1048576.; // 2^20

const int nRepeat=1000;

typedef quantumoperator::Tridiagonal<RANK> Tridiagonal ;
typedef quantumoperator::Tridiagonal<   1> Tridiagonal1;

typedef quantumdata::StateVector<RANK> StateVector;

typedef Tridiagonal1::Diagonal Diagonal1;


int main()
{
  int dim=mathutils::round(pow(totalDim,1./RANK));
  
  Diagonal1 zero(dim), minus(dim-2), plus(dim-2);

  {
    using randomized::fillWithRandom;
    fillWithRandom(plus,fillWithRandom(minus,fillWithRandom(zero)));
  }

  const Tridiagonal1 element(zero,2,minus,plus);


#define PRINT_factor(z,n,data) element BOOST_PP_IIF(BOOST_PP_LESS(n, BOOST_PP_SUB(RANK, 1) ) , *, ) 
  
  const Tridiagonal tridiag( BOOST_PP_REPEAT(RANK,PRINT_factor,~) );

#undef  PRINT_factor

  cerr<<"Dimensions: "<<element.getTotalDimension()<<' '<<tridiag.getTotalDimension()//<<' '<<minus.size()<<' '<<plus.size()
      <<endl;


  StateVector psi(tridiag.getDimensions()), dpsidt(psi);

  {
    using randomized::fillWithRandom;
    fillWithRandom(psi(),1002);
  }

  {
    boost::progress_timer t;

    for (int i=0; i<nRepeat; i++)
      tridiag.apply(psi(),dpsidt());
  }

}
