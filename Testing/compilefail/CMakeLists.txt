# Copyright Raimar Sandner 2012–2020.
# Copyright András Vukics 2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

add_compilefail_test(fail_Composite "Composite not consistent" CPPQED::CPPQEDcore)

add_compilefail_test(fail_DensityOperator "Incorrect number of subscripts for DensityOperator" CPPQED::CPPQEDcore)

add_compilefail_test(fail_StateVector "Incorrect number of subscripts for StateVector" CPPQED::CPPQEDcore)
