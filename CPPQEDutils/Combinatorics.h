// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDCORE_UTILS_COMBINATORICS_H_INCLUDED
#define   CPPQEDCORE_UTILS_COMBINATORICS_H_INCLUDED

#include "Conversions.h"
#include "MathExtensions.h"

#include <array>
#include <numeric>
#include <stdexcept>
#include <vector>


namespace cppqedutils {


template<size_t NoOfSites>
class CWR_Dir // Combinations with repetitions
{
public:
  using Configuration=std::array<size_t,NoOfSites>;
  
  using Impl=std::vector<Configuration>;

  CWR_Dir(size_t k) : impl_{[=] () {
    Impl res{longDouble2Size_t(multiChoose<long double>(NoOfSites,k))};
    
    auto recurse=[&](typename Impl::iterator i, size_t actualSite, size_t remainingParticles, auto& recurseFunction) -> void {
      if (actualSite==NoOfSites-1) { (*i)[actualSite]=remainingParticles;}
      else {
        for (size_t putHere=0; putHere<=remainingParticles; ++putHere) {
          recurseFunction(i,actualSite+1,remainingParticles-putHere,recurseFunction);
          auto offset=longDouble2Ptrdiff_t(multiChoose<long double>(NoOfSites-actualSite-1,remainingParticles-putHere));
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
  // NoOfSites is the number of objects from which you can choose and k is the number to be chosen
  // E.g. for bosonic states on a lattice: NoOfSites is the number of lattice sites, k is the number of bosonic particles.

  const Configuration& operator[](size_t i) const {return impl_[i];}
  
  struct SubscriptingException : public std::range_error
  {
    SubscriptingException(const Configuration& c) : std::range_error("Configuration not found"), conf(c) {}
    
    const Configuration conf;
  };

  size_t operator[](const Configuration& c) const
  {
    if (auto res=std::find(impl_.begin(),impl_.end(),c); res==impl_.end() ) throw SubscriptingException{c};
    else return res-impl_.begin();
  }
  
  const Impl& operator()() const {return impl_;}

private:
  const Impl impl_;

};



} // cppqedutils


#endif // CPPQEDCORE_UTILS_COMBINATORICS_H_INCLUDED
