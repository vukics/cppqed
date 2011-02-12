// -*- C++ -*-
#ifndef _INDEXER_IMPLEMENTATIONS_H
#define _INDEXER_IMPLEMENTATIONS_H

#include <boost/fusion/sequence/intrinsic/at_c.hpp>

#include<boost/preprocessor/iteration/iterate.hpp>
#include<boost/preprocessor/repetition.hpp>

#define BOOST_PP_ITERATION_LIMITS (1,11)
#define BOOST_PP_FILENAME_1 "details/IndexerImplementationsSpecialization.h"

#include BOOST_PP_ITERATE()

#undef BOOST_PP_FILENAME_1
#undef BOOST_PP_ITERATION_LIMITS

#endif // _INDEXER_IMPLEMENTATIONS_H
