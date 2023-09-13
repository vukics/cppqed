#include <array>
#include <cstddef>
#include <iostream>


constexpr size_t calculateRANK(int)
{
  return 1;
}


template <size_t RANK>
constexpr size_t calculateRANK(std::array<int,RANK>)
{
  return RANK;
}


template <typename T> constexpr bool isOffsetType(T) {return false;}

constexpr bool isOffsetType(int) {return true;}

template <size_t RANK> constexpr bool isOffsetType(std::array<int,RANK>) {return true;}


template <auto... offsets> requires (isOffsetType(offsets) && ...)
struct Tridiagonal
{
  static constexpr size_t RANK=calculateRANK(std::array{offsets...}.front());
};



int main()
{
  Tridiagonal<1,2,1,3> a;
  Tridiagonal<std::array{1,2},std::array{1,3}> b;
  std::cerr<<decltype(a)::RANK<<" "<<decltype(b)::RANK<<std::endl;
}
