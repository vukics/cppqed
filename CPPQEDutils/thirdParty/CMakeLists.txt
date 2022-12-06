# Copyright András Vukics 2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

add_library(ThirdParty INTERFACE)

target_include_directories(ThirdParty INTERFACE Xoshiro-cpp pcg-cpp/include popl/include)
