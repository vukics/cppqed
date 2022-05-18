#include "MultiArray.h"

#include <iostream>
#include <random>
#include <tuple>

using namespace cppqedutils;
using namespace boost::json;

int main()
{
  std::vector<double> vec(1000);
  
  std::span<double> span(vec);
  
  std::cout<<vec.size()<<" "<<span.size()<<std::endl;
  
  MultiArray<double,5> ma{{10,11,13,9,8}};
  
  double& v1=ma(1,5,6,3,4) ;
  
  std::cout<<v1<<std::endl;
  
  // std::cout<<ma(10,5,6,3,4)<<std::endl;
  // double& v2=ma(1,5,6,std::tuple<int>{3},4) ;

  // cppqedutils::MultiArray<const double,5> mac; // "std::vector must have a non-const, non-volatile value_type"

  // cppqedutils::MultiArrayView<const double,5> mavc;

  // const double& v3=mavc(1,5,6,3,4) ;
  // double& v4=mavc(1,5,6,3,4) ;
  
  // double& v5=ma(1,5,6,3) ;
  // double& v6=ma(1,5,6,3,4,2) ;
  
  auto offsets(multiarray::calculateSlicesOffsets(ma.extents,ma.strides,tmptools::vector<2,4,0>));
  
  for (auto o : offsets) std::cout<<o<<std::endl;
  
  auto sr=slicesRange(ma,tmptools::vector<2,4,0>,offsets);

  std::uniform_real_distribution distro{1.,10.};
  std::mt19937 gen{1001};
  
  MultiArray<double,3> maSmall{{5,7,4},[&] (size_t s) {
    MultiArray<double,3>::StorageType data(s);
    for (auto& d : data) d=distro(gen);
    return data;
  }};  
  
  std::string maSmallSerialized=serialize( value_from( maSmall ) );
  
  std::cout << maSmallSerialized << std::endl;
  
  MultiArray<double,3> maSmallReconstructed{ value_to<MultiArray<double,3>>(parse(maSmallSerialized)) };
  
  std::cout << (maSmall==maSmallReconstructed) << std::endl;
  
  for (auto o : maSmallReconstructed.dataView) std::cout<<o<<" "/*<<std::endl*/;
  
  std::cout<<std::endl;
  
  std::cout<<serialize( value_from( maSmallReconstructed ) )<<std::endl;
  
  std::cout<< (maSmallSerialized == serialize( value_from( maSmallReconstructed ) ) ) << std::endl;
  
//   {
//     MultiArray<double,5> ma{{2,3,4,2,4}};
//     
//     for (size_t i=0; i<ma.dataView.size(); ++i) ma.dataView[i]=i;
//     
//     for (auto&& slice : slicesRange(ma,tmptools::vector<2,4,0>,
//                                     multiarray::calculateSlicesOffsets(ma.extents,ma.strides,tmptools::vector<2,4,0>)))
//       std::cout<<slice(0,0,0)<<std::endl;
//   }
  
}
