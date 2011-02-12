///////////////////////////////////////////////////////////////////////////////
/// \file filter.hpp
///   Defines the filter range adaptor, as well as the make_filter_range() helper
//
///////////////////////////////////////////////////////////////////////////////
#ifndef RANGE_FILTER_EN_12_09_2004_HPP
#define RANGE_FILTER_EN_12_09_2004_HPP

#include <boost/static_assert.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include <boost/range/result_iterator.hpp>
#include "./iterator_cast.hpp"
#include "./detail/adaptor_base.hpp"
#include "./detail/adl_begin_end.hpp"

namespace boost
{

/// \brief filters a range using the specified predicate
///
/// filters a range using the specified predicate
///
template<typename FwdRng,typename Pred>
struct filter_range
  : boost::iterator_range<
        boost::filter_iterator<
            Pred
          , BOOST_DEDUCED_TYPENAME boost::range_result_iterator<FwdRng>::type
        >
    >
{
    typedef BOOST_DEDUCED_TYPENAME boost::range_result_iterator<FwdRng>::type   base_iterator;
    typedef boost::filter_iterator<Pred,base_iterator>                          iterator;
    typedef boost::iterator_range<iterator>                                     base_range;

    explicit filter_range(FwdRng &rng, Pred pred)
      : base_range(
            boost::make_filter_iterator(pred,range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng))
          , boost::make_filter_iterator(pred,range_ex_detail::adl_end(rng),range_ex_detail::adl_end(rng))
        )
    {
    }
};

/// \brief filters a range using the specified predicate
///
/// filters a range using the specified predicate
///
template<typename FwdRng,typename Pred>
filter_range<FwdRng,Pred> make_filter_range(FwdRng& rng,Pred pred)
{
    return filter_range<FwdRng,Pred>(rng,pred);
}

/// \overload
template<typename FwdRng,typename Pred>
filter_range<FwdRng const,Pred> make_filter_range(FwdRng const& rng,Pred pred)
{
    return filter_range<FwdRng const,Pred>(rng,pred);
}


/// filter_range_adaptor
///   INTERNAL ONLY
struct filter_range_adaptor
{
    template<typename Rng,typename Args>
    struct apply
    {
        BOOST_STATIC_ASSERT((boost::tuples::length<Args>::value==1));
        typedef filter_range<
            Rng
          , BOOST_DEDUCED_TYPENAME boost::tuples::element<0,Args>::type
        > type;
    };

    template<typename Rng,typename Args>
    static BOOST_DEDUCED_TYPENAME apply<Rng,Args>::type
    make_range(Rng & rng,Args const & args)
    {
        return make_filter_range(rng, boost::tuples::get<0>(args));
    }
};

namespace adaptor
{
    namespace
    {
        /// \brief the filter range adaptor
        ///
        /// the filter range adaptor
        ///
        range_adaptor<filter_range_adaptor> const filter;
    }
}

} // namespace boost

#endif
