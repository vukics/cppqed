#include "Sigma.h"

#include "Mode.h"

using namespace std;
using namespace mode;

using blitz::shape;

using quantumdata::Sigma;

int main()
{
  Sigma<4,1> sigma4_1;
  Sigma<2,3> sigma2_3;

  Tridiagonal a1(aop(10)), a2(aop(8).dagger());
  
  structure::Types<4>::StateVectorLow psi(shape(6,10,8,5)), dpsidt(shape(6,10,8,5));

  psi=0; dpsidt=0; psi(1,7,6,3)=1;

  (sigma4_1*(a1*a2)*sigma2_3).apply<4>(psi,dpsidt);

  cout<<dpsidt(4,6,7,2)<<endl;

  /*
  {

    (sigma4_1*sigma4_1*sigma4_1).apply<3>(psi,dpsidt);
  }


  {
    structure::Types<4>::StateVectorLow psi, dpsidt;

    (sigma4_1*(a*a)*sigma4_1).apply<4>(psi,dpsidt);

  }
  */
  return 0;

}
