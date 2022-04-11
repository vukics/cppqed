// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines i/oarchive types depending on the #BZ_HAVE_BOOST_SERIALIZATION macro}
#ifndef   CPPQEDCORE_UTILS_ARCHIVE_H_INCLUDED
#define   CPPQEDCORE_UTILS_ARCHIVE_H_INCLUDED

#ifndef BZ_HAVE_BOOST_SERIALIZATION

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
/** \see \refBoost{Boost.Serialization,serialization/doc/archive_reference.html#trivial} */
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
/** \copydetails trivial_iarchive */
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

#else // BZ_HAVE_BOOST_SERIALIZATION

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

namespace cppqedutils {

typedef boost::archive::binary_iarchive iarchive; ///< delegated to \refBoost{Boost.Serialization,serialization}
typedef boost::archive::binary_oarchive oarchive; ///< \copydoc iarchive

} // cppqedutils

#endif // BZ_HAVE_BOOST_SERIALIZATION

#endif // CPPQEDCORE_UTILS_ARCHIVE_H_INCLUDED
