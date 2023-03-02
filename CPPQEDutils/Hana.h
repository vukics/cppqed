// Copyright András Vukics 2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include <boost/hana.hpp>

namespace hana=boost::hana;

template <typename S> concept hana_sequence = hana::Sequence<S>::value;
