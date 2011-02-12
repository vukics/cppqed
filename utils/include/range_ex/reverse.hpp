///////////////////////////////////////////////////////////////////////////////
/// \file reverse.hpp
///   Defines the reverse range adaptor, as well as the make_reverse_range() helper
//
///////////////////////////////////////////////////////////////////////////////
#ifndef RANGE_REVERSE_EN_12_09_2004_HPP
#define RANGE_REVERSE_EN_12_09_2004_HPP

#include <boost/static_assert.hpp>
#include <boost/iterator/reverse_iterator.hpp>
#include <boost/range/result_iterator.hpp>
#include "./iterator_cast.hpp"
#include "./detail/adaptor_base.hpp"
#include "./detail/adl_begin_end.hpp"

namespace boost
{

/// \brief a range that is the reverse of the original
///
/// a range that is the reverse of the original
///
template<typename BidiRng>
struct reverse_range
  : boost::iterator_range<
        boost::reverse_iterator<
            BOOST_DEDUCED_TYPENAME boost::range_result_iterator<BidiRng>::type
        >
    >
{
    typedef BOOST_DEDUCED_TYPENAME boost::range_result_iterator<BidiRng>::type  base_iterator;
    typedef boost::reverse_iterator<base_iterator>                              iterator;
    typedef boost::iterator_range<iterator>                                     base_range;

    explicit reverse_range(BidiRng &rng)
      : base_range(
            boost::make_reverse_iterator(range_ex_detail::adl_end(rng))
          , boost::make_reverse_iterator(range_ex_detail::adl_begin(rng))
        )
    {
    }
};

/// \brief produces a range that is the reverse of the original
///
/// produces a range that is the reverse of the original
///
template<typename BidiRng>
reverse_range<BidiRng> make_reverse_range(BidiRng &rng)
{
    return reverse_range<BidiRng>(rng);
}

/// \overload
template<typename BidiRng>
reverse_range<BidiRng const> make_reverse_range(BidiRng const &rng)
{
    return reverse_range<BidiRng const>(rng);
}

/// reverse_range_adaptor
///   INTERNAL ONLY
struct reverse_range_adaptor
{
    template<typename Rng,typename Args>
    struct apply
    {
        BOOST_STATIC_ASSERT((boost::tuples::length<Args>::value==0));
        typedef reverse_range<Rng> type;
    };

    template<typename Rng,typename Args>
    static BOOST_DEDUCED_TYPENAME apply<Rng,Args>::type
    make_range(Rng & rng, Args)
    {
        return make_reverse_range(rng);
    }
};

namespace adaptor
{
    namespace
    {
        /// \brief the reverse range adaptor
        ///
        /// the reverse range adaptor
        ///
        range_adaptor<reverse_range_adaptor> const reverse;
    }
}

} // namespace boost

#endif
