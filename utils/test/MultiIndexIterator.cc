#include "MultiIndexIterator.tcc"

#include <boost/lambda/lambda.hpp>
namespace bll=boost::lambda;
namespace mpl=boost::mpl;


#include <algorithm>

using namespace std;
using namespace cpputils;

int main()
{
  {
    typedef MultiIndexIterator<10> MII10;
    typedef MII10::IdxTiny IdxTiny;

    const IdxTiny 
      lbound(0,0,0,0,0,0,0,0,0,0),
      ubound(23,4,9,8,32,7,4,9,3,28);

    MII10 begin(lbound,ubound,mpl::bool_<false>()), end(lbound,ubound,mpl::bool_<true>());

    // cout<<*(++begin)<<' '<<*end<<endl;
  }

  {
    typedef MultiIndexIterator<4> MII4;
    typedef MII4::IdxTiny IdxTiny;

    const IdxTiny 
      lbound(0,0,0,0),
      ubound(5,4,9,8);

    MII4 begin(lbound,ubound,mpl::bool_<false>()), end(lbound,ubound,mpl::bool_<true>());

    // cout<<*(++begin)<<' '<<*end<<endl;

    for_each(begin,end,cout<<bll::_1<<' ');

  }

  
  // for_each
}
