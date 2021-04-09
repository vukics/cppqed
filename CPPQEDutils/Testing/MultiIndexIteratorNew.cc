// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "MultiIndexIterator.h"

using namespace cppqedutils;

namespace mpl=boost::mpl;

typedef IdxTiny<3> Idx;



int main()
{
  Idx lbound(-2,0,1), ubound(-1,1,3);

  MultiIndexIterator<3> 
    begin(lbound,ubound,mii::begin),
    end(lbound,ubound,mii::end);

  return !(all(*begin++==lbound) && 
           all(*begin++==Idx(-2,0,2)) &&
           all(*begin++==Idx(-2,0,3)) &&
           all(*begin++==Idx(-2,1,1)) &&
           all(*begin++==Idx(-2,1,2)) &&
           all(*begin++==Idx(-2,1,3)) &&
           all(*begin++==Idx(-1,0,1)) &&
           all(*begin++==Idx(-1,0,2)) &&
           all(*begin++==Idx(-1,0,3)) &&
           all(*begin++==Idx(-1,1,1)) &&
           all(*begin++==Idx(-1,1,2)) &&
           all(*begin++==Idx(-1,1,3)) &&
           begin++==end && all(*end==Idx(0,0,1)) && all(*begin==Idx(0,0,2)));

}
