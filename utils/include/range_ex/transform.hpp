///////////////////////////////////////////////////////////////////////////////
/// \file transform.hpp
///   Defines the transform range adaptor, as well as the make_transform_range() helper
//
///////////////////////////////////////////////////////////////////////////////
#ifndef RANGE_TRANSFORM_EN_12_09_2004_HPP
#define RANGE_TRANSFORM_EN_12_09_2004_HPP

#include <boost/static_assert.hpp>
#include <boost/iterator/transform_iterator.hpp>
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
template<typename FwdRng,typename UnaryFunc>
struct transform_range
  : boost::iterator_range<
        boost::transform_iterator<
            UnaryFunc
          , BOOST_DEDUCED_TYPENAME boost::range_result_iterator<FwdRng>::type
        >
    >
{
    typedef BOOST_DEDUCED_TYPENAME boost::range_result_iterator<FwdRng>::type   base_iterator;
    typedef boost::transform_iterator<UnaryFunc,base_iterator>                  iterator;
    typedef boost::iterator_range<iterator>                                     base_range;

    explicit transform_range(FwdRng &rng, UnaryFunc fun)
      : base_range(
            boost::make_transform_iterator(range_ex_detail::adl_begin(rng), fun)
          , boost::make_transform_iterator(range_ex_detail::adl_end(rng), fun)
        )
    {
    }
};

/// \brief transforms a range using the specified unary function
///
/// transforms a range using the specified unary function
///
template<typename FwdRng,typename UnaryFunc>
transform_range<FwdRng,UnaryFunc> make_transform_range(FwdRng& rng, UnaryFunc fun)
{
    return transform_range<FwdRng,UnaryFunc>(rng, fun);
}

/// \overload
template<typename FwdRng,typename UnaryFunc>
transform_range<FwdRng const,UnaryFunc> make_transform_range(FwdRng const& rng, UnaryFunc fun)
{
    return transform_range<FwdRng const,UnaryFunc>(rng, fun);
}

/// transform_range_adaptor
///   INTERNAL ONLY
struct transform_range_adaptor
{
    template<typename Rng,typename Args>
    struct apply
    {
        BOOST_STATIC_ASSERT((boost::tuples::length<Args>::value==1));
        typedef transform_range<
            Rng
          , BOOST_DEDUCED_TYPENAME boost::tuples::element<0,Args>::type
        > type;
    };

    template<typename Rng,typename Args>
    static BOOST_DEDUCED_TYPENAME apply<Rng,Args>::type
    make_range(Rng & rng, Args const & args)
    {
        return make_transform_range(rng, boost::tuples::get<0>(args));
    }
};

///////////////////////////////////////////////////////////////////////////////
// transform
//
namespace adaptor
{
    namespace
    {
        /// \brief the transform range adaptor
        ///
        /// the transform range adaptor
        ///
        range_adaptor<transform_range_adaptor> const transform;
    }
}

} // namespace boost

#endif
