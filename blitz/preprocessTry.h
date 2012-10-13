#include "array-impl.h"

#define DEFAULT_print(z, n, data) nilArraySection()
#define BOOST_PP_ITERATION_LIMITS (1,BLITZ_ARRAY_LARGEST_RANK)
#define BOOST_PP_FILENAME_1 "blitz/details/SlicingOperatorReentrant.h"

#include BOOST_PP_ITERATE()

#undef BOOST_PP_FILENAME_1
#undef BOOST_PP_ITERATION_LIMITS
#undef DEFAULT_print
