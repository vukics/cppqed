#include <iostream>

#include "MultiArrayComplex.h"

#include "Random.h"

#include <tuple>

using namespace cppqedutils;
using namespace boost::json;

constexpr auto ra20741 = retainedAxes<2,0,7,4,1>;


int main()
{
  MultiArray<dcomp,9> array{{5,2,3,2,4,3,2,4,6}};

  Extents<5> 
    filteredExtents{multiarray::filterIn<ra20741>(array.extents)}/*,
    filteredStrides{multiarray::filterIn<ra20741>(array.strides)}*/;

  std::cerr<<array.dataView.size()<<std::endl;

  {
    std::uniform_real_distribution d{-1.0,1.0};
    std::mt19937_64 re{1001};
    std::ranges::generate(array.mutableView().dataView, [&]() {return std::complex{d(re),d(re)};});
  }

  std::cerr<<array(1,1,2,2,0,1,1,3,2)<<std::endl<<serialize( value_from( filteredExtents ) )<<std::endl;
  
  {
    std::cerr<<"*** Simple range: ***"<<std::endl;
    
    auto range{sliceRangeSimple<ra20741>(array)};
    
    auto iter{range.begin()};
    
    std::cerr<<serialize( value_from( iter->extents ) )<<std::endl;
    for (size_t i=0; i<10; (++i, ++iter));
    
    std::cerr<<(*iter)(2,2,1,3,0)<<std::endl;
    
    auto rangeRecurse{sliceRangeSimple<retainedAxes<1,3>>(*iter)};

    std::cerr<<serialize( value_from( rangeRecurse.begin()->extents ) )<<std::endl;
    
    auto iter2{rangeRecurse.begin()};
    for (size_t i=0; i<3; (++i, ++iter2));

    std::cerr<<(*iter2)(3,2)<<std::endl;
    
  }
  
  // auto offsets{multiarray::calculateSlicesOffsets<ra20741>(array.extents,array.strides)};
  
  // for (SliceIterator<std::vector<size_t>,dcomp,5> iter{filteredExtents,filteredStrides,offsets.begin(),array.dataView}; iter!=offsets.end(); ++iter)
  //  std::cerr<<(*iter)(1,0,2,1,0)<<std::endl;
  
  // SliceRangeReferencing<dcomp,5> sliceRange{filteredExtents,filteredStrides,offsets,array.dataView};
  
  // for (MultiArrayView<dcomp,5> mav : sliceRange) std::cerr<<mav(1,0,2,1,0)<<std::endl;
  
  // for (MultiArrayView<dcomp,5> macv : sliceRange) std::cerr<<macv.dataView[0]/*(1,0,2,1,0)*/<<std::endl;
  
  {
    std::cerr<<"*** Proxy range owning: ***"<<std::endl;
    
    auto range{sliceRange<ra20741>(array)};
    
    auto iter{range.begin()};
    
    std::cerr<<serialize( value_from( iter->extents ) )<<std::endl;
    for (size_t i=0; i<10; (++i, ++iter));
    
    std::cerr<<(*iter)(2,2,1,3,0)<<std::endl;
    
    auto rangeRecurse{sliceRange<retainedAxes<1,3>>(*iter)};

    std::cerr<<serialize( value_from( rangeRecurse.begin()->extents ) )<<std::endl;

    auto iter2{rangeRecurse.begin()};
    for (size_t i=0; i<3; (++i, ++iter2));

    std::cerr<<(*iter2)(3,2)<<std::endl;

  }

  {
    std::cerr<<"*** Proxy range referencing: ***"<<std::endl;
    
    auto offsets{calculateSlicesOffsets<ra20741>(array.extents)};
    
    auto range{sliceRange<ra20741>(array,offsets)};
    
    auto iter{range.begin()};
    
    std::cerr<<serialize( value_from( iter->extents ) )<<std::endl;
    for (size_t i=0; i<10; (++i, ++iter));
    
    std::cerr<<(*iter)(2,2,1,3,0)<<std::endl;
    
    auto rangeRecurse{sliceRange<retainedAxes<1,3>>(*iter)};

    std::cerr<<serialize( value_from( rangeRecurse.begin()->extents ) )<<std::endl;
    
    auto iter2{rangeRecurse.begin()};
    for (size_t i=0; i<3; (++i, ++iter2));

    std::cerr<<(*iter2)(3,2)<<std::endl;

  }

  {
    MultiArray<dcomp,4> a1{{5,2,3,2}};
    MultiArray<dcomp,5> a2{{4,3,2,4,6}};
  
    {
      std::uniform_real_distribution d{-1.0,1.0};
      std::mt19937_64 re{1001};
      std::ranges::generate(a1.mutableView().dataView, [&]() {return std::complex{d(re),d(re)};} );
      std::ranges::generate(a2.mutableView().dataView, [&]() {return std::complex{d(re),d(re)};} );
    }

    auto res{directProduct<4,5>(a1,a2)};

  }
  
}



void f()
{
  std::vector<double> vec(1000);
  
  std::span<double> span(vec);
  
  std::cout<<vec.size()<<" "<<span.size()<<std::endl;
  
  MultiArray<double,5> ma{{10,11,13,9,8}};

  {
    MultiArrayConstView<double,5> macv{ma};
    // macv(2,3,1,6,2)=1.1;
  }
  
  std::cout<<ma.dataView.size()<<std::endl;
  
  double& v1=ma(1,5,6,3,4) ;
  
  std::cout<<v1<<std::endl;
  
  // std::cout<<ma(10,5,6,3,4)<<std::endl;
  // double& v2=ma(1,5,6,std::tuple<int>{3},4) ;

  // cppqedutils::MultiArray<const double,5> mac; // "std::vector must have a non-const, non-volatile value_type"

  cppqedutils::MultiArrayView<const double,5> mavc;

  // const double& v3=mavc(1,5,6,3,4) ;
  // double& v4=mavc(1,5,6,3,4) ;
  
  // double& v5=ma(1,5,6,3) ;
  // double& v6=ma(1,5,6,3,4,2) ;
  
  auto offsets(multiarray::calculateSlicesOffsets<retainedAxes<2,4,0>>(ma.extents,ma.strides));
  
  for (auto o : offsets) std::cout<<o<<std::endl;
  
  auto sr=sliceRange<retainedAxes<2,4,0>>(ma,offsets);

  for (auto mav : sr) std::cout<<mav.offset<<std::endl;
  
  MultiArray<double,3> maSmall{{5,7,4},[&] (size_t s) {
    MultiArray<double,3>::StorageType data(s);
    std::uniform_real_distribution uniform{1.,10.};
    std::mt19937 gen{1001};

    for (auto& d : data) d=uniform(gen);
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
  
  {
    Extents<5> idx{0,0,0,0,0}, ext{5,2,4,3,2};
    for (size_t i=0; i<multiarray::calculateExtent(ext); (++i, incrementMultiIndex(idx, ext) ) ) std::cout<<serialize(value_from(idx))<<std::endl;
  }
  
  std::cout<<hana::Sequence<decltype(compileTimeRange<0,10>)>::value<<hana::Sequence<decltype(retainedAxes<1,5,10>)>::value<<std::endl;
  
  std::cout<<serialize(value_from(std::make_tuple(1.2,"hello",32ul)))<<std::endl;
  
//   {
//     MultiArray<double,5> ma{{2,3,4,2,4}};
//     
//     for (size_t i=0; i<ma.dataView.size(); ++i) ma.dataView[i]=i;
//     
//     for (auto&& slice : slicesRange(ma,retainedAxes<2,4,0>,
//                                     multiarray::calculateSlicesOffsets(ma.extents,ma.strides,retainedAxes<2,4,0>)))
//       std::cout<<slice(0,0,0)<<std::endl;
//   }

  // std::views::join({a1.extents,a2.extents})
  
  {
    std::vector<size_t> v1{1,2,3,4}, v2{5,6,7};
    auto joined(std::views::join(std::vector{v1,v2}));
    
    std::cerr<<serialize( value_from( std::vector<size_t>(joined.begin(),joined.end()) ) )<<std::endl;
  }
  
}



auto calcIndex=[]<size_t RANK> (Extents<RANK> extents, auto... i) requires (sizeof...(i)==RANK) {
  size_t res=0;
  auto e=extents.begin();

  ( [&] {
#ifndef   NDEBUG
    if (i >= *e) throw std::range_error("Index position: "+std::to_string(e-extents.begin())+", index value: "+std::to_string(i)+", extent: "+std::to_string(*e));
#endif // NDEBUG
    return res+=(*e++)*i;} (), ... );
  
  return res;
  
};


void g()
{
  MultiArray<dcomp,9> array{{5,2,3,2,4,3,2,4,6}};

  Extents<5> 
    filteredExtents{multiarray::filterIn<ra20741>(array.extents)},
    filteredStrides{multiarray::filterIn<ra20741>(array.strides)};

  std::cerr<<calcIndex(filteredExtents,1,2,3,2,0)<<std::endl;
}
