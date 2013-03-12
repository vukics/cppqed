#ifndef   UTILS_INCLUDE_ARCHIVE_H_INCLUDED
#define   UTILS_INCLUDE_ARCHIVE_H_INCLUDED

#ifdef    DO_NOT_USE_BOOST_SERIALIZATION

#include <cstddef> // std::size_t
#include <iosfwd>

#include <boost/mpl/bool.hpp>


namespace cpputils {

//////////////////////////////////////////////////////////////
// class trivial_(i/o)archive 
// cf. http://www.boost.org/doc/libs/1_53_0/libs/serialization/doc/archive_reference.html#trivial

class trivial_iarchive {
public:
  typedef boost::mpl::bool_<false> is_saving; 
  typedef boost::mpl::bool_<true> is_loading;
  
  trivial_iarchive(std::istream&) {}
  
  template<class T> void register_type(){}
  template<class T> trivial_iarchive & operator>>(const T & t){
      return *this;
  }
  template<class T> trivial_iarchive & operator&(const T & t){
      return *this >> t;
  }
  void load_binary(void *address, std::size_t count){};
};


class trivial_oarchive {
public:
  typedef boost::mpl::bool_<true> is_saving; 
  typedef boost::mpl::bool_<false> is_loading;

  trivial_oarchive(std::ostream&) {}

  template<class T> void register_type(){}
  template<class T> trivial_oarchive & operator<<(const T & t){
      return *this;
  }
  template<class T> trivial_oarchive & operator&(const T & t){
      return *this << t;
  }
  void save_binary(void *address, std::size_t count){};
};


typedef trivial_iarchive iarchive;
typedef trivial_oarchive oarchive;

} // cpputils

#else  // DO_NOT_USE_BOOST_SERIALIZATION

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

namespace cpputils {

typedef boost::archive::binary_iarchive iarchive;
typedef boost::archive::binary_oarchive oarchive;

} // cpputils

#endif // DO_NOT_USE_BOOST_SERIALIZATION

#endif // UTILS_INCLUDE_ARCHIVE_H_INCLUDED
