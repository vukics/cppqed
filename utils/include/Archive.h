#ifndef   UTILS_INCLUDE_ARCHIVE_H_INCLUDED
#define   UTILS_INCLUDE_ARCHIVE_H_INCLUDED

#ifdef    DO_NOT_USE_BOOST_SERIALIZATION

#include <iosfwd>

namespace cpputils {

typedef std::istream iarchive;
typedef std::ostream oarchive;

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
