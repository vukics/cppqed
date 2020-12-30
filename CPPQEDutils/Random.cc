// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

#include "Random.h"


using namespace std;

// Everything that is *defined* in GSL must be declared here, since GSL is linked only into CPPQEDutils

cpputils::GSL_RandomEngine::GSL_RandomEngine(result_type s/*, const gsl_rng_type* ran_gen_type*/) 
  : ranGen_(gsl_rng_alloc(gsl_rng_taus2),[](auto p) {gsl_rng_free(p);}),
    name_(std::string(gsl_rng_name(ranGen_.get())))
{
  gsl_rng_set(ranGen_.get(),s);
}


void cpputils::GSL_RandomEngine::seed(result_type value) {gsl_rng_set(ranGen_.get(),value);}

auto cpputils::GSL_RandomEngine::operator()() -> result_type {return gsl_rng_get(ranGen_.get());}

auto cpputils::GSL_RandomEngine::min() const -> result_type {return gsl_rng_min(ranGen_.get());}

auto cpputils::GSL_RandomEngine::max() const -> result_type {return gsl_rng_max(ranGen_.get());}

void cpputils::GSL_RandomEngine::write(std::ostream& os) const
{
  (os/*<<name_*/).write(static_cast<const char*>(gsl_rng_state(ranGen_.get())),gsl_rng_size(ranGen_.get()));
}


void cpputils::GSL_RandomEngine::read(istream& is)
{
  /* string n; is>>n;
  if (n!=name_) {
    is.clear(ios_base::failbit);
    throw std::runtime_error("GSL_RandomEngine archive load -- wrong implementation ID, expected "+name_+", found "+n);
  }*/
  is.read(static_cast<char*>(gsl_rng_state(ranGen_.get())),gsl_rng_size(ranGen_.get()));
}
