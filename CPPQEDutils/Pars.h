// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "Traits.h"

#include "popl.hpp"

inline std::string parsedCommandLine="";

popl::OptionParser optionParser(std::string pre="", std::string post="");

void parse(popl::OptionParser&, int argc, const char* const argv[]);

template <typename T>
auto& add(popl::OptionParser& op, std::string option, std::string description, T defaultValue, auto&& binding)
{
  op.add<popl::Value<T>>("",option,description,defaultValue,std::forward<decltype(binding)>(binding));
  return op;
}

template <typename T>
auto& add(popl::OptionParser& op, std::string option, std::string mod, std::string description, T defaultValue, auto&& binding)
{
  return add(op,option+mod,description,defaultValue,std::forward<decltype(binding)>(binding));
}


inline auto& add(popl::OptionParser& op, std::string option, std::string description, bool* binding)
{
  op.add<popl::Switch>("",option,description,binding);
  return op;
}


inline auto& add(popl::OptionParser& op, std::string option, std::string mod, std::string description, bool* binding)
{
  return add(op,option+mod,description,binding);
}



/// TODO: insert this more directly into popl
inline auto& addTitle(popl::OptionParser& op, std::string title, std::string mod = "")
{
  op.add<popl::Value<std::string>>("","\n### "+title+mod,"###","###");
  return op;
}


namespace parameters {


template <typename T, std::convertible_to<T> U>
auto _(std::string option, std::string description, const U& defaultValue, T& binding)
{
  return std::tuple<std::string,std::string,T,T&>{option,description,defaultValue,binding};
}


/// `T` is either a string (title), or a tuple
template <typename... T> requires ( ... && ( decltype( ::cppqedutils::multilambda {
  [] <typename U, typename V> (const std::tuple<std::string,std::string,U,V&> &) { return std::is_convertible<U,V>{}; },
  [] <typename U> (const U &) { return std::is_convertible<U,std::string>{}; }
  } (std::declval<T>()))::value ) )
popl::OptionParser& add_dispatch(std::string mod, popl::OptionParser& op, const T&... t)
{
  ::cppqedutils::multilambda worker {
    [&op,mod] <typename U, typename V> (const std::tuple<std::string,std::string,U,V&> & t ) {
      ::add(op,std::get<0>(t),mod,std::get<1>(t),std::get<2>(t),&std::get<3>(t));
    },
    [&op,mod] (const std::string& title) {::addTitle(op,title,mod);}
  };

  (worker(t), ...);

  return op;
}


} // parameters



popl::OptionParser& addTuple(std::string mod, popl::OptionParser& op, const auto&... t) {return parameters::add_dispatch(mod,op,t...);}

popl::OptionParser& addTuple(popl::OptionParser& op, const auto&... t) {return parameters::add_dispatch("",op,t...);}



// popl::OptionParser& addTuple(popl::OptionParser& op, const auto&... t)
// {
//   return addTuple(op,"",t...);
// }
