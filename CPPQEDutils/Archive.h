// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/json.hpp>

namespace cppqedutils {

typedef boost::archive::binary_iarchive iarchive; ///< delegated to \refBoost{Boost.Serialization,serialization}
typedef boost::archive::binary_oarchive oarchive; ///< \copydoc iarchive


std::string toStringJSON(auto&& v) {return ::boost::json::serialize( ::boost::json::value_from( std::forward<decltype(v)>(v) ) ) ;}


// template <typename T, typename Archive> concept serializable = requires (Archive& ar, T& t) {serialize(ar,t)};
// the problem with this is that the Archive type has to be explicitly specified when using the concept
// cf. Q&As here https://stackoverflow.com/a/63452199/1171157 and https://stackoverflow.com/a/67092860/1171157


} // cppqedutils
