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

const size_t nRepeat=10000000;


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
  // size_t res=0;
  switch(rank_) {
  case 1:
    return 0;
  case 2:
    return 0;
  case 3:
    return 0;
  case 4:
    return 0;
  case 5:
#ifndef NDEBUG
    if (idx[0]>=sizes_[0] || idx[1]>=sizes_[1] || idx[2]>=sizes_[2] || idx[3]>=sizes_[3] || idx[4]>=sizes_[4]) abort();
#endif
    return idx[0]*strides_[0]+idx[1]*strides_[1]+idx[2]*strides_[2]+idx[3]*strides_[3]+idx[4]*strides_[4];
  case 6:
    return 0;
  case 7:
    return 0;
  case 8:
    return 0;
  case 9:
    return 0;
  case 10:
#ifndef NDEBUG
    if (idx[0]>=sizes_[0] || idx[1]>=sizes_[1] || idx[2]>=sizes_[2] || idx[3]>=sizes_[3] || idx[4]>=sizes_[4] || idx[5]>=sizes_[5] || idx[6]>=sizes_[6] || idx[7]>=sizes_[7] || idx[8]>=sizes_[8] || idx[9]>=sizes_[9]) abort();
#endif
    return 
      idx[0]*strides_[0]+idx[1]*strides_[1]+idx[2]*strides_[2]+idx[3]*strides_[3]+idx[4]*strides_[4]+
      idx[5]*strides_[5]+idx[6]*strides_[6]+idx[7]*strides_[7]+idx[8]*strides_[8]+idx[9]*strides_[9];
  default :
    abort();
  }

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

  {

    Size size; size+=16,15,14,17,18;
    Data temp(16*15*14*17*18);
    Size idx(5), idxv(5);
    /*

    { // Handcrafted
      RunTimeRanked target(size,fillRan(temp));
      const RunTimeRanked source(size,fillRan(temp));

      {
	boost::progress_timer t;
	for (size_t count=0; count<nRepeat; ++count)
	  for (int j=0; j<16; ++j) {
	    idx[0]=13; idx[1]=4; idx[2]=11; idx[3]=j; idx[4]=12;
	    idxv[0]=11; idxv[1]=6; idxv[2]=10; idxv[3]=j+1; idxv[4]=11;
	    target(idx,source(idxv));
	  }
      }
    }
    { // Blitz++
      Array<double,5> target(16,15,14,17,18), source(16,15,14,17,18);
      fillRan(target); fillRan(source);

      {
	boost::progress_timer t;
	for (size_t count=0; count<nRepeat; ++count)
	  for (int j=0; j<16; ++j)
	    target(13,4,11,j,12)=source(11,6,10,j+1,11);
      }
    }
    */
  }


  double  arrayBI[1028160
    ];
  double* arrayDY=new double[1028160
    ];


  {
    boost::progress_timer t;
    for (size_t count=0; count<nRepeat; ++count)
      for (int j=0; j<601243
	     ; ++j) {
	arrayBI[123405+j]=arrayBI[296732+j];
	arrayDY[296732+j]=arrayBI[123405+j];
      }
  }

  {
    boost::progress_timer t;
    for (size_t count=0; count<nRepeat; ++count)
      for (int j=0; j<601243
	     ; ++j)
	arrayDY[123405+j]=arrayDY[296732+j];
  }


}
