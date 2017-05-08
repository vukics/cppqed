// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
/// \briefFile{Collecting the macros and for compilation configuration to one place for documentation with doxygen + documents directories. This file is never included in the framework.}

/// \dir utils \brief General purpose utilities, which can be considered a small but rather diverse library in themselves

/// \dir composites \brief Comprises composite systems BinarySystem and Composite and relating utilities

#define DO_NOT_USE_BOOST_SERIALIZATION 0 ///< Governs whether the \refBoost{Boost.Serialization,serialization} library is installed and can/should be used

#define DO_CONSIDER_EXPLICITLY_SPECIALIZED_TRIDIAGONAL_APPLIES 1 ///< Governs which implementation of quantumoperator::Tridiagonal::apply is considered.
