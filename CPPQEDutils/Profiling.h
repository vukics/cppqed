// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

#undef  PROGRESS_TIMER_OUT_POINT
#undef  PROGRESS_TIMER_IN_POINT

#ifdef    USE_BOOST_PROGRESS_TIMER_PROFILING

#include <boost/progress.hpp>

#define PROGRESS_TIMER_IN_POINT(OS) { std::ostream& progressReportStream=OS; boost::progress_timer boost_progress_timer(progressReportStream);

#define PROGRESS_TIMER_OUT_POINT(REPORT)  progressReportStream<<REPORT<<": "; }

#else

#define PROGRESS_TIMER_IN_POINT(OS) {

#define PROGRESS_TIMER_OUT_POINT(REPORT) }

#endif // USE_BOOST_PROGRESS_TIMER_PROFILING


#undef    USE_BOOST_PROGRESS_TIMER_PROFILING
