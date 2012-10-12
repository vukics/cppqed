#include "impl/ComplexArrayExtensions.tcc"

#include <boost/preprocessor/iteration/iterate.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/arithmetic/mul.hpp>
#include <boost/preprocessor/arithmetic/div.hpp>

#define BOOST_PP_ITERATION_LIMITS (1,BOOST_PP_DIV(BLITZ_ARRAY_LARGEST_RANK,2))
#define BOOST_PP_FILENAME_1 "../src/details/HCH_ImplementationsSpecialization.h"

#include BOOST_PP_ITERATE()

#undef BOOST_PP_FILENAME_1
#undef BOOST_PP_ITERATION_LIMITS

