# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

add_library(CPPQEDcore_quantumdata SHARED DistributionFunctions.cc LazyDensityOperatorFFT.cc) # Transformation.cc off for the moment
add_library(CPPQED::CPPQEDcore_quantumdata ALIAS CPPQEDcore_quantumdata)

target_link_libraries(CPPQEDcore_quantumdata PUBLIC CPPQED::CPPQEDutils)
