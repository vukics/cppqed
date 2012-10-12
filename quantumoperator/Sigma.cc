#include "Sigma.h"

#include "TMP_Tools.h"

using quantumdata::Types;
using blitz::Range;

#include<boost/preprocessor/iteration/iterate.hpp>
#include<boost/preprocessor/repetition/enum.hpp>
#include<boost/preprocessor/arithmetic/sub.hpp>

#define BOOST_PP_ITERATION_LIMITS (2,BOOST_PP_SUB(BLITZ_ARRAY_LARGEST_RANK,1))
#define BOOST_PP_FILENAME_1 "../quantumoperator/details/PartialProjectImplementationSpecializations.h"

#include BOOST_PP_ITERATE()

#undef BOOST_PP_FILENAME_1
#undef BOOST_PP_ITERATION_LIMITS

