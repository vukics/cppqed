# Copyright András Vukics 2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

add_library(ThirdPartyRandom INTERFACE)

target_include_directories(ThirdPartyRandom INTERFACE Xoshiro-cpp pcg-cpp/include popl/include)
