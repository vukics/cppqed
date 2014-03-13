// -*- C++ -*-
/// \briefFile{Documents directories. This file is never included in the framework.}

/// \dir CPPQEDscripts \brief Comprises the supported scripts of the framework. Supported means that they will always remain an integral part of the framework’s distribution and are tested by the testsuite(s).

#define DO_NOT_USE_BOOST_SERIALIZATION 0 ///< Governs whether the \refBoost{Boost.Serialization,serialization} library is installed and can/should be used

#define DO_CONSIDER_EXPLICITLY_SPECIALIZED_TRIDIAGONAL_APPLIES 1 ///< Governs which implementation of quantumoperator::Tridiagonal::apply is considered.
