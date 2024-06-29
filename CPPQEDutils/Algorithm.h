// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Range algorithms not yet in STL in C++20}
#pragma once

#include <algorithm>
#include <iterator>
#include <numeric>
#include <ranges>
#include <vector>

#include <boost/math/special_functions/binomial.hpp>

namespace cppqedutils {


/// std::array`s of different size cannot be std::views::join`ed because the size is part of the type
/** this solution based on C++17 fold expressions comes from [here](http://stackoverflow.com/a/42774523/1171157) */
template <typename Type, std::size_t... sizes>
constexpr auto concatenate(const std::array<Type, sizes>&... arrays)
{
  std::array<Type, (sizes + ...)> result;
  std::size_t index{};

  ((std::copy_n(arrays.begin(), sizes, result.begin() + index), index += sizes), ...);

  return result;
}



template <class T>
inline T multiChoose(unsigned n, unsigned k) {return boost::math::binomial_coefficient<T>(n+k-1,k);}


// Combinations with repetitions
// NoOfSites is the number of objects from which you can choose and k is the number to be chosen
// E.g. for bosonic states on a lattice: NoOfSites is the number of lattice sites, k is the number of bosonic particles.
template<size_t NoOfSites>
struct CWR_Dir
{
public:
  using Configuration=std::array<size_t,NoOfSites>;

  using Impl=std::vector<Configuration>;

  CWR_Dir(size_t k) : configurations{[=] {
    Impl res{static_cast<size_t>(multiChoose<double>(NoOfSites,k))};

    auto recurse=[&](typename Impl::iterator i, size_t actualSite, size_t remainingParticles, auto& recurseFunction) -> void {
      if (actualSite==NoOfSites-1) { (*i)[actualSite]=remainingParticles;}
      else {
        for (size_t putHere=0; putHere<=remainingParticles; ++putHere) {
          recurseFunction(i,actualSite+1,remainingParticles-putHere,recurseFunction);
          auto offset=static_cast<ptrdiff_t>(multiChoose<double>(NoOfSites-actualSite-1,remainingParticles-putHere));
          for (auto j=i; j<i+offset; ++j) (*j)[actualSite]=putHere;
          i+=offset;
        }
      }
    };

    recurse(res.begin(),0,k,recurse);

#ifndef   NDEBUG
    for (const auto& config : res) if (std::accumulate(config.begin(),config.end(),0) != k) throw std::logic_error("Problem in CWR configurations");
#endif // NDEBUG

    return res;

  }() } {}

  const Configuration& operator[](size_t i) const {return configurations[i];}

  struct SubscriptingException : public std::range_error
  {
    SubscriptingException(const Configuration& c) : std::range_error("Configuration not found"), conf(c) {}

    const Configuration conf;
  };

  size_t operator[](const Configuration& c) const
  {
    if (auto res=std::find(configurations.begin(),configurations.end(),c); res==configurations.end() ) throw SubscriptingException{c};
    else return res-configurations.begin();
  }

  size_t operator()(std::convertible_to<size_t> auto ... i) const requires (sizeof...(i)==NoOfSites) {return operator[]({i...});}

  const Impl configurations;

};


} // cppqedutils
