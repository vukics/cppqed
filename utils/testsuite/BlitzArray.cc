#include "BlitzArray.h"

#include "Randomized.h"
// #include<memcopy>

#include "Range.h"

#include<boost/bind.hpp>

using namespace std;
using namespace randomized;

typedef blitz::Array<dcomp,6> Array;

typedef blitzplusplus::Array<dcomp,6> MyArray;


int main()
{
  Randomized::SmartPtr Ran(MakerGSL()(1001));

  Array array(blitz::shape(2,4,3,2,5,4));

  boost::generate(array,bind(&Randomized::dcompRan,Ran));

  array.transposeSelf(1,3,0,4,2,5);

  Array arrayv(array.copy()), arrayvv(array.shape()); arrayvv=array; 
  
  cout<<array.ordering()<<arrayv.ordering()<<arrayvv.ordering()<<endl;

  cout<<all(arrayvv==array)<<endl;

  dcomp* data=new dcomp[array.size()];

  memcpy(data,array.dataFirst(),array.size());

  blitzplusplus::Array<dcomp,6> myArray(MyArray(array,blitzplusplus::ShallowCopy()).clone(data),blitzplusplus::ShallowCopy());

  cout<<all((2.*myArray)!=array)<<endl;

  // return 0;

}

