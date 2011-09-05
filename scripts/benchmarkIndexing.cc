#include "Randomized.h"
#include "Range.h"

#include <blitz/array.h>

#include <boost/assign/list_of.hpp>
#include <boost/assign/std/vector.hpp>

#include <boost/progress.hpp>

#include <vector>


using namespace std;
using namespace randomized;
using namespace boost::assign;

using blitz::Array;

typedef vector<double> Data;
typedef vector<int> Size;

const size_t nRepeat=1000000;


class RunTimeRanked
{
public:
  static const Size calculateStrides(const Size& sizes)
  {
    Size res(1,1);
    for (Size::const_iterator i=sizes.begin(); i!=sizes.end()-1; ++i)
      res.push_back((*i)*(*--res.end()));
    return res;
  }

  RunTimeRanked(const Size& sizes, const Data& data)
    : rank_(sizes.size()),
      data_(data),
      sizes_(sizes),
      strides_(calculateStrides(sizes))
  {
    // for (Size::const_iterator i=strides_.begin(); i!=strides_.end(); ++i) cout<<*i<<' '<<endl;
  }

  double operator()(const Size&) const;
  
  void operator()(const Size&, double);

private:
  size_t calculateIndex(const Size&) const;

  const size_t rank_;

  Data data_;

  const Size sizes_, strides_;

};


double RunTimeRanked::operator()(const Size& idx) const
{
  return data_[calculateIndex(idx)];
}


void RunTimeRanked::operator()(const Size& idx, double v)
{
  data_[calculateIndex(idx)]=v;
}


size_t RunTimeRanked::calculateIndex(const Size& idx) const
{
  size_t res=0;
  for (Size::const_iterator
#ifndef NDEBUG
	 z=sizes_.begin(),
#endif
	 i=idx.begin(),
	 s=strides_.begin();
       i!=idx.end();
#ifndef NDEBUG
       ++z,
#endif
	 ++i, ++s) {
#ifndef NDEBUG
    if (!(*z>*i)) abort();
#endif
    res+=(*i)*(*s);
  }
  return res;
}


template<typename A>
const A& fillRan(A& data)
{
  Randomized::SmartPtr Ran(MakerGSL()(1001));
  boost::generate(data,bind(&Randomized::operator(),Ran));
  return data;
}



int main()
{
  /////////
  // RANK=5
  /////////

  { // Handcrafted
    Data temp(16*15*14*17*18);

    RunTimeRanked target(list_of(16)(15)(14)(17)(18),fillRan(temp));

    const RunTimeRanked source(list_of(16)(15)(14)(17)(18),fillRan(temp));

    {
      boost::progress_timer t;
      for (size_t count=0; count<nRepeat; ++count)
	for (int j=0; j<16; ++j)
	  target(list_of(13)(4)(11)(j)(12),source(list_of(11)(6)(10)(j+1)(11)));
    }
  }
  { // Blitz++
    Array<double,5> target(16,15,14,17,18), source(16,15,14,17,18);
    fillRan(target); fillRan(source);

    {
      boost::progress_timer t;
      for (size_t count=0; count<100*nRepeat; ++count)
	for (int j=0; j<16; ++j)
	  target(13,4,11,j,12)=source(11,6,10,j+1,11);
    }
  }


}
