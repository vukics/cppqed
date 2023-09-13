// #include <boost/json.hpp>

#include "MultiArray.h"

#include <boost/hana.hpp>

#include <complex>
#include <concepts>
#include <iostream>
#include <ranges>
#include <set>
#include <vector>


/*
struct CustomIterator {
  using wrapped_type = std::vector<int>::iterator;

  using value_type = wrapped_type::value_type;
  using difference_type = wrapped_type::difference_type;
  using pointer = wrapped_type::pointer;
  using reference = wrapped_type::reference;
  using iterator_category = std::input_iterator_tag;

  wrapped_type iter;

  reference operator*() const { return *iter; }
  pointer operator->() const { return &*iter; }

  CustomIterator& operator++() { ++iter; return *this; }

  friend auto operator<=>(const CustomIterator& lhs, const CustomIterator& rhs) = default;

};


struct CustomRange {
  using wrapped_type = std::vector<int>;

  wrapped_type vector;

  auto begin() const { return vector.begin(); }

  auto end() const { return vector.end(); }
};
*/


struct CustomIterator {
  using value_type = size_t;
  using difference_type = ptrdiff_t;
/*  using pointer = size_t;   // ???
  using reference = size_t; // ???*/
  using iterator_category = std::input_iterator_tag;

  mutable value_type idx;

  size_t operator*() const { return idx; }

  CustomIterator& operator++() { ++idx; return *this; }

  friend auto operator<=>(const CustomIterator& lhs, const CustomIterator& rhs) = default;

  // endpoint guard
  friend auto operator!=(const CustomIterator& lhs, size_t rhs) {return lhs.idx!=rhs;}
};


struct CustomRange {
  size_t extent;

  CustomIterator begin() const { return CustomIterator{0}; }

  size_t end() const { return extent; }
};


static_assert( std::ranges::range<CustomRange> );
// static_assert( std::ranges::input_range<CustomRange> );


namespace hana = boost::hana;


constexpr auto ra{cppqedutils::retainedAxes<2,0,7,4,1>};

constexpr auto sorted([s(ra)] mutable { std::ranges::sort(s); return s; } () );

constexpr auto uniquedSorted([us(sorted)] mutable { std::ranges::unique(us); return us; } () );

constexpr bool hasDuplicate=[us(ra)] mutable { std::ranges::sort(us); const auto [first, last] = std::ranges::unique(us); return (first==last); } () ;

static_assert( hasDuplicate );

// static_assert( std::ranges::equal(uniquedSorted,sorted) );

//static_assert( std::ranges::equal(std::ranges::unique(sorted), sorted ) );


struct A
{
  A(const A&) = delete;
  A(A&&) {std::cerr<<"move ctor called"<<std::endl;}
  A() {std::cerr<<"default ctor called"<<std::endl;}
};

struct B {};
struct C {};


constexpr auto typeOptions = hana::tuple_t<A,B>;


template <typename T>
concept type_options = std::same_as<T,A> || std::same_as<T,B>;


template <typename L>
concept functional = hana::fold( typeOptions, true, [] <typename E> (bool s, E) { return s &&  ( requires (const L& l, typename E::type rho) { { l(rho) } -> std::convertible_to<double> ; } ) ; } );



auto f = [] (type_options auto) {return 1.;};

// static_assert(functional_helper<decltype(f),A>);
// static_assert(functional<decltype(f)>);


A g() {return A();}



template <size_t RANK>
struct CompareOffsets
{
  using Offsets = std::array<int,RANK>;

  bool operator() (const Offsets& o1, const Offsets o2) const {
    auto recurse = [&] (Offsets::const_iterator i1, Offsets::const_iterator i2, const auto r) {
      if (i1==o1.end()) return false;
      else if (*i1==*i2) return r(++i1,++i2,r);
      else return *i1<*i2;
    };
    return recurse(o1.begin(),o2.begin(),recurse);
  }
};



int main()
{
  std::vector v{2,4,6,8,10}, vv{3,6,9,12,15};

  CustomRange r{10};

  for (auto i : r) std::cerr<<i<<std::endl;

//  if (auto iter=std::ranges::find(r,7uz); iter!=std::end(r)) std::cerr<<*iter<<std::endl;


 // for (auto&& [i,j] : std::views::zip(CustomRange{10},CustomRange{10}) ) std::cerr<<i<<" "<<j<<std::endl;

/*  boost::json::object o, oo;

  o["int"]=2;
  o["double"]=2.1;//boost::json::value_from(std::complex{1.2,-3.1});
  o["vector"]=boost::json::value_from(std::vector{1,3,2});

  oo["object"]=o;
  oo["pi"]=3.14;

  std::cerr<<oo<<std::endl;*/
/*
  A a(g()), aa(std::move(g()));

  std::set<CompareOffsets<5>::Offsets,CompareOffsets<5>> s{
    {-2,3,7,2,-1},
    {-3,1,2,-4,5},
    {-2,3,4,1,-5},
    {1,3,2,1,0},
    {1,3,2,1,2},
    {0,3,-2,-6,1},
    {-3,1,2,-4,5}
  };

  for (const auto& o : s) std::cerr<<cppqedutils::json(o)<<std::endl;*/

}
