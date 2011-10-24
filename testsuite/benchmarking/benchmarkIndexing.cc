#include "Randomized.h"
#include "Range.h"

#include <blitz/array.h>

#include <boost/assign/std/vector.hpp>

#include <boost/progress.hpp>

#include <vector>


using namespace boost::assign;

using blitz::Array;

typedef std::vector<double> Data;
typedef std::vector<int   > Size;

typedef blitz::TinyVector<int,5> Idx;

const size_t nRepeat=100000000;


class RunTimeRanked // In this form works only for Rank=5
{
public:
  static const Size calculateStrides(const Size& sizes)
  {
    Size res(1,1);
    for (Size::const_iterator i=sizes.begin(); i!=sizes.end()-1; ++i) res+=(*i)*(*--res.end());
    return res;
  }

  RunTimeRanked(const Size& sizes, const Data& data)
    : rank_(sizes.size()), data_(data), sizes_(sizes), strides_(calculateStrides(sizes)), stridesRaw_(&strides_[0]) {}


  double operator()(int i0, int i1, int i2, int i3, int i4) const     {return data_[calculateIndex(i0,i1,i2,i3,i4)]  ;}
  void   operator()(int i0, int i1, int i2, int i3, int i4, double v) {       data_[calculateIndex(i0,i1,i2,i3,i4)]=v;}
  // ... similar brute-force implementations for all ranks.

private:
  size_t calculateIndex(int i0, int i1, int i2, int i3, int i4) const {return i0*strides_[0]+i1*strides_[1]+i2*strides_[2]+i3*strides_[3]+i4*strides_[4];}
  // ... similar brute-force implementations for all ranks.

  const size_t rank_;
  Data data_;
  const Size sizes_, strides_;

  const Idx stridesRaw_;

};


template<typename A>
const A& fillRan(A& data)
{
  using namespace randomized;
  Randomized::SmartPtr Ran(MakerGSL()(1001));
  boost::generate(data,bind(&Randomized::operator(),Ran));
  return data;
}



int main()
{
  /////////
  // RANK=5
  /////////

  double d=0;

  { // Handcrafted
    Size size; size+=16,15,14,17,18 ;
    Data temp       (16*15*14*17*18);
    RunTimeRanked       target(size,fillRan(temp));
    const RunTimeRanked source(size,fillRan(temp));

    {
      boost::progress_timer t;
      for (size_t count=0; count<nRepeat; ++count)
	for (int j=1; j<16; ++j) {
	  // target(13,4,11,j,12,source(11,6,10,j+1,11));
	  d+=target(13,4,11,j-1,12);
	}
    }
  }
  { // Blitz++
    Array<double,5> target(16,15,14,17,18), source(16,15,14,17,18);
    fillRan(target); fillRan(source);
    
    {
      boost::progress_timer t;
      for (size_t count=0; count<nRepeat; ++count)
	for (int j=1; j<16; ++j) {
	  // target(13,4,11,j,12)=source(11,6,10,j+1,11);
	  d+=target(13,4,11,j-1,12);
	}
      
    }
  }
  std::cout<<d<<std::endl;

}
