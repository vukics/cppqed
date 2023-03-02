#include <boost/hana.hpp>

#include <iostream>

namespace hana = boost::hana;

// this is only for consistency
template <typename S> concept hana_sequence = hana::Sequence<S>::value;

template <typename T>
concept callable_with_int = requires(T&& t, int i) { t(i); };

/*
// spurious helper needed?
template <typename T>
constexpr bool sequenceOfCallablesHelper = false;

template <typename... T>
constexpr bool sequenceOfCallablesHelper<hana::tuple<T...>> =  hana::fold(hana::tuple_t<T...>, true, [] (bool s, auto element) {
    return s && callable_with_int<typename decltype(element)::type>;
  });
*/

template <typename S>
concept sequence_of_callables = ( hana_sequence<S> && !!hana::all_of(
    decltype(hana::transform(std::declval<S>(), hana::typeid_)){},
    []<class T>(T) { return callable_with_int<typename T::type>; })
/* sequenceOfCallablesHelper<S>*/ );
// why the !! ?
// why the decltype around hana::transform ?

void applyCallable(const sequence_of_callables auto& ha, int i)
{
  hana::for_each(ha, [&] (auto e) {std::cerr<<"inside apply – "; e(i);});
}

int main()
{
  applyCallable(hana::make_tuple(
    [] (int i) {std::cerr<<"inside lambda0 "<<2*i<<"!\n";},
    [] (int i) {std::cerr<<"inside lambda1 "<<i<<"!\n";}/*,
    [] () {std::cerr<<"inside lambda2!\n";}*/ // this would spoil compilation, as expected
  ),4);
}



namespace {

template <size_t> struct A {};


template <typename H, size_t RANK>
concept callable_with_A = requires(H&& h, A<RANK> a) { h(a); };

template <typename H, size_t RANK>
concept callable_with_A_and_int = requires(H&& h, A<RANK> b, int i) { { h(b,i) } -> std::convertible_to<int>; };

template <typename H, size_t RANK>
concept callable = ( callable_with_A<H,RANK> || callable_with_A_and_int<H,RANK> ) ;


template <size_t RANK, typename T>
constexpr bool sequenceOfCallablesHelper = false;

template <size_t RANK, typename... T>
constexpr bool sequenceOfCallablesHelper<RANK,hana::tuple<T...>> =  hana::fold(hana::tuple_t<T...>, true, [] (bool s, auto element) {
    return s && callable<typename decltype(element)::type,RANK>;
  });


template <typename S, size_t RANK>
concept sequence_of_callables = ( hana_sequence<S> && sequenceOfCallablesHelper<RANK,S> );


template <size_t RANK>
auto applyCallable(const sequence_of_callables<RANK> auto& ha, A<RANK> a, int i)
{
  hana::for_each(ha, [&] (auto e) {
    if constexpr (callable_with_A<decltype(e),RANK>) {std::cerr<<"CWA "; e(a);}
    else if (callable_with_A_and_int<decltype(e),RANK>) {std::cerr<<"CWAI "; e(a,i);}
  });
}


int foo()
{
  applyCallable(hana::make_tuple(
    [] (A<2>, int i) {std::cerr<<"inside "<<i<<"!\n"; return i;},
    [] (A<2>) {std::cerr<<"inside!\n";},
    [] (A<2>) {std::cerr<<"inside other!\n";}
  ),A<2>{},4);
}

}

#include <optional>
#include <vector>

void g()
{
  auto ht{hana::make_tuple("string",1)}; static_assert(hana_sequence<decltype(ht)>);

  std::vector v{1,2,3}; static_assert(!hana_sequence<decltype(v)>);
}

#include "MultiArrayComplex.h"

#include "Random.h"

#include <tuple>

/*
using namespace cppqedutils;
using namespace boost::json;
*/

void f()
{
  auto v=boost::json::value_from(std::array{1,2,3});
  
  cppqedutils::Extents<3> /*std::array*/ a{1,2,3};
  std::cerr<<cppqedutils::toStringJSON(a)<<std::endl;

  cppqedutils::MultiArray<int,9> array{{2,4,6}};
  std::cerr<<cppqedutils::toStringJSON(array)<<std::endl;

  
}
