///////////////////////////////////////////////////////////////////////////////
/// \file indirect.hpp
///   Defines the indirect range adaptor, as well as the make_indirect_range() helper
//
///////////////////////////////////////////////////////////////////////////////
#ifndef RANGE_INDIRECT_EN_12_09_2004_HPP
#define RANGE_INDIRECT_EN_12_09_2004_HPP

#include <boost/static_assert.hpp>
#include <boost/iterator/indirect_iterator.hpp>
#include <boost/range/result_iterator.hpp>
#include "./iterator_cast.hpp"
#include "./detail/adaptor_base.hpp"
#include "./detail/adl_begin_end.hpp"

namespace boost
{

/// \brief produced an indirect range
///
/// produced an indirect range, where each element in the new
/// range is the result of a dereference of the corresponding
/// element in the original range.
///
template<typename FwdRng>
struct indirect_range
  : boost::iterator_range<
        boost::indirect_iterator<
            BOOST_DEDUCED_TYPENAME boost::range_result_iterator<FwdRng>::type
        >
    >
{
    typedef BOOST_DEDUCED_TYPENAME boost::range_result_iterator<FwdRng>::type   base_iterator;
    typedef boost::indirect_iterator<base_iterator>                             iterator;
    typedef boost::iterator_range<iterator>                                     base_range;

    explicit indirect_range(FwdRng &rng)
      : base_range(
            boost::make_indirect_iterator(range_ex_detail::adl_begin(rng))
          , boost::make_indirect_iterator(range_ex_detail::adl_end(rng))
        )
    {
    }
};

/// \brief produced an indirect range
///
/// produced an indirect range, where each element in the new
/// range is the result of a dereference of the corresponding
/// element in the original range.
///
template<typename FwdRng>
indirect_range<FwdRng> make_indirect_range(FwdRng& rng)
{
    return indirect_range<FwdRng>(rng);
}

/// \overload
template<typename FwdRng>
indirect_range<FwdRng const> make_indirect_range(FwdRng const& rng)
{
    return indirect_range<FwdRng const>(rng);
}

/// indirect_range_adaptor
///   INTERNAL ONLY
struct indirect_range_adaptor
{
    template<typename Rng,typename Args>
    struct apply
    {
        BOOST_STATIC_ASSERT((boost::tuples::length<Args>::value==0));
        typedef indirect_range<Rng> type;
    };

    template<typename Rng,typename Args>
    static BOOST_DEDUCED_TYPENAME apply<Rng,Args>::type
    make_range(Rng & rng, Args)
    {
        return make_indirect_range(rng);
    }
};

namespace adaptor
{
    namespace
    {
        /// \brief the indirect range adaptor
        ///
        /// the indirect range adaptor
        ///
        range_adaptor<indirect_range_adaptor> const indirect;
    }
}

} // namespace boost

#endif
