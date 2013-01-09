// -*- C++ -*-
#ifndef   UTILS_INCLUDE_FACTORY_H_INCLUDED
#define   UTILS_INCLUDE_FACTORY_H_INCLUDED


// NEEDS_WORK: generalize this for an arbitrary number of parameters with Boost.Preprocessor

#define CPPUTILS_FACTORY_TEMPLATE_PARAMETERS typename CFTP1, typename CFTP2, typename CFTP3
#define CPPUTILS_FACTORY_CONSTRUCTOR_PARAMETERS_WITH_DEFAULTS const CFTP1& cftp1=cpputils::factory::_1(), const CFTP2& cftp2=cpputils::factory::_2(), const CFTP3& cftp3=cpputils::factory::_3()
#define CPPUTILS_FACTORY_CONSTRUCTOR_PARAMETERS const CFTP1& cftp1, const CFTP2& cftp2, const CFTP3& cftp3
#define CPPUTILS_FACTORY_CONSTRUCTOR_ARGUMENTS cftp1,cftp2,cftp3

namespace cpputils {


namespace factory {

// Placeholders:
class _1 {};
class _2 {};
class _3 {};
  
} // factory
  
template<typename Base>
class Factory : public Base
{
public:
  template<typename A1, typename A2, typename A3>
  Factory(const A1& a1, const A2& a2, const A3& a3) : Base(a1,a2,a3) {}
  
  template<typename A1, typename A2>
  Factory(const A1& a1, const A2& a2, factory::_3) : Base(a1,a2) {}
  
  template<typename A1>
  Factory(const A1& a1, factory::_2, factory::_3) : Base(a1) {}

  Factory(factory::_1, factory::_2, factory::_3) : Base() {}
  
};


} // cpputils


#endif // UTILS_INCLUDE_FACTORY_H_INCLUDED
