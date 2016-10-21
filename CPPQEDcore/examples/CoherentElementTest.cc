// Copyright András Vukics 2006–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include <MathExtensions.h>
#include <Pars.tcc>

#include <boost/math/special_functions/factorials.hpp>

using namespace std;
using namespace boost::math;
using namespace mathutils;

int main(int argc, char* argv[])
{
  parameters::ParameterTable p;
  dcomp         & alpha=p.add("alpha","",dcomp(-1,2));
  unsigned long & max  =p.add("max"  ,"",200ul      );
  
  update(p,argc,argv,"--");
  
  cerr<<max_factorial<double>::value<<" "<<max_factorial<long double>::value<<endl;
  for (unsigned long n=1; n<max; ++n) {
    const dcomp
      straight(n<max_factorial<double>::value ? pow(alpha,n)/sqrt(factorial<double>(n)) : dcomp()),
      stirling(n ? pow(2*n*PI,-.25)*pow(alpha/sqrt(n/EULER),n) : dcomp(1.));
    cout<<n<<" "<<straight.real()<<" "<<straight.imag()<<" "<<stirling.real()<<" "<<stirling.imag()<<" "<<abs((straight-stirling)/(straight+stirling))<<endl;
  }
}