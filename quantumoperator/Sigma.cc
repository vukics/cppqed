#include "Sigma.h"

using quantumdata::Types;
using blitz::Range;

#include<boost/preprocessor/iteration/iterate.hpp>
#include<boost/preprocessor/repetition.hpp>
#include<boost/preprocessor/arithmetic/sub.hpp>

#define BOOST_PP_ITERATION_LIMITS (2,10)
#define BOOST_PP_FILENAME_1 "../quantumoperator/details/PartialProjectImplementationSpecializations.h"

#include BOOST_PP_ITERATE()

#undef BOOST_PP_FILENAME_1
#undef BOOST_PP_ITERATION_LIMITS

