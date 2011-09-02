#include "Randomized.h"

#include <blitz/array.h>

#include <boost/assign/list_of.hpp>
#include <boost/assign/std/vector.hpp>

#include <boost/progress.hpp>

#include <vector>


using namespace std;
using namespace randomized;
using namespace boost::assign;


typedef vector<double> Data;
typedef vector<size_t> Size;


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
    for (Size::const_iterator i=strides_.begin(); i!=strides_.end(); ++i) cout<<*i<<' '<<endl;
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
	 i=idx.begin(),
	 s=strides_.begin();
       i!=idx.end();
       ++i, ++s)
    res+=(*i)*(*s);
  return res;
}


int main()
{
  /////////
  // RANK=5
  /////////

  {
    Data temp(16*15*14*17*18);

    {
      Randomized::SmartPtr Ran(MakerGSL()(1001));
      boost::generate(temp,bind(&Randomized::dcompRan,Ran));
    }

    RunTimeRanked rtr(list_of(16)(15)(14)(17)(18),temp);
  }
}
