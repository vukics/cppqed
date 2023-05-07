// #include <boost/json.hpp>

#include <boost/hana.hpp>

#include <complex>
#include <iostream>
#include <vector>


namespace hana = boost::hana;


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


int main()
{
/*  boost::json::object o, oo;

  o["int"]=2;
  o["double"]=2.1;//boost::json::value_from(std::complex{1.2,-3.1});
  o["vector"]=boost::json::value_from(std::vector{1,3,2});

  oo["object"]=o;
  oo["pi"]=3.14;

  std::cerr<<oo<<std::endl;*/

  A a(g()), aa(std::move(g()));

}
