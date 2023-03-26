#include <boost/json.hpp>

#include <complex>
#include <iostream>
#include <vector>

int main()
{
  boost::json::object o, oo;

  o["int"]=2;
  o["double"]=2.1;//boost::json::value_from(std::complex{1.2,-3.1});
  o["vector"]=boost::json::value_from(std::vector{1,3,2});

  oo["object"]=o;
  oo["pi"]=3.14;

  std::cerr<<oo<<std::endl;

}
