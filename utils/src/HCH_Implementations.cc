#include "ComplexArrayExtensions.h"

#include<boost/preprocessor/iteration/iterate.hpp>
#include<boost/preprocessor/repetition.hpp>
#include<boost/preprocessor/arithmetic/mul.hpp>

#define BOOST_PP_ITERATION_LIMITS (1,5)
#define BOOST_PP_FILENAME_1 "details/HCH_ImplementationsSpecialization.h"

#include BOOST_PP_ITERATE()

#undef BOOST_PP_FILENAME_1
#undef BOOST_PP_ITERATION_LIMITS
