// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

namespace cppqedutils {

typedef boost::archive::binary_iarchive iarchive; ///< delegated to \refBoost{Boost.Serialization,serialization}
typedef boost::archive::binary_oarchive oarchive; ///< \copydoc iarchive


// template <typename T, typename Archive> concept serializable = requires (Archive& ar, T& t) {serialize(ar,t)};
// the problem with this is that the Archive type has to be explicitly specified when using the concept
// cf. Q&As here https://stackoverflow.com/a/63452199/1171157 and https://stackoverflow.com/a/67092860/1171157


} // cppqedutils


/*

#include <cstddef> // size_t
#include <iosfwd>


namespace cppqedutils {

//////////////////////////////////////////////////////////////

class trivial_archive_base {
public:
  typedef std::false_type is_saving ;
  typedef std:: true_type is_loading;

  template<class T> void register_type(){}
  
};

/// Trivial iarchive disabling serialization
// \see \refBoost{Boost.Serialization,serialization/doc/archive_reference.html#trivial}
class trivial_iarchive : public trivial_archive_base {
public:
  
  trivial_iarchive(std::istream&) {}
  
  template<class T> trivial_iarchive & operator>>(const T & ){
      return *this;
  }
  template<class T> trivial_iarchive & operator&(const T & t){
      return *this >> t;
  }
  void load_binary(void *, size_t){};
};


/// Trivial oarchive disabling serialization
class trivial_oarchive : public trivial_archive_base {
public:

  trivial_oarchive(std::ostream&) {}

  template<class T> trivial_oarchive & operator<<(const T & ){
      return *this;
  }
  template<class T> trivial_oarchive & operator&(const T & t){
      return *this << t;
  }
  void save_binary(void *, size_t){};
};


typedef trivial_iarchive iarchive; ///< delegated to trivial_iarchive
typedef trivial_oarchive oarchive; ///< delegated to trivial_oarchive

} // cppqedutils

*/
