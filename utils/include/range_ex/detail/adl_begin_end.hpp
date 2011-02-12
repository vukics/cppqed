///////////////////////////////////////////////////////////////////////////////
//
// adl_begin_end.hpp
//
/////////////////////////////////////////////////////////////////////////////

#if defined(_MSC_VER) && _MSC_VER >= 1000
#pragma once
#endif

#ifndef ADL_BEGIN_END_EN_14_12_2004
#define ADL_BEGIN_END_EN_14_12_2004

#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/result_iterator.hpp>

namespace boost
{
    namespace range_ex_detail
    {
        template<typename Rng>
        inline BOOST_DEDUCED_TYPENAME boost::range_result_iterator<Rng>::type
        adl_begin(Rng & rng)
        {
	  return boost::begin(rng);
        }

        template<typename Rng>
        inline BOOST_DEDUCED_TYPENAME boost::range_result_iterator<Rng>::type
        adl_end(Rng & rng)
        {
	  return boost::end(rng);
        }
    }
}

#endif
