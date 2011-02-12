///////////////////////////////////////////////////////////////////////////////
//
// File: adaptor_base.hpp
//
///////////////////////////////////////////////////////////////////////////////
#ifndef ADAPTOR_BASE_EN_12_08_2004_HPP
#define ADAPTOR_BASE_EN_12_08_2004_HPP

#include <boost/iterator/iterator_traits.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/tuple/tuple.hpp>

namespace boost
{

///////////////////////////////////////////////////////////////////////////////
// range_adaptor
//
template<typename Base, typename Args = boost::tuple<> >
struct range_adaptor
    : Base
{
    range_adaptor()
        : Base()
        , args()
    {
    }

    range_adaptor(Base const & base, Args const & args = Args())
        : Base(base)
        , args(args)
    {
    }

    range_adaptor(range_adaptor<Base,Args> const &that)
        : Base(that)
        , args(that.args)
    {
    }

    range_adaptor<Base,Args> const & operator ()() const
    {
        return *this;
    }

    template<typename Arg1>
    range_adaptor<Base,boost::tuple<Arg1> > operator ()(Arg1 const & arg1) const
    {
        return range_adaptor<Base,boost::tuple<Arg1> >(*this,boost::make_tuple(arg1));
    }

    template<typename Arg1,typename Arg2>
    range_adaptor<Base,boost::tuple<Arg1,Arg2> > operator ()(Arg1 const & arg1,Arg2 const & arg2) const
    {
        return range_adaptor<Base,boost::tuple<Arg1,Arg2> >(*this,boost::make_tuple(arg1,arg2));
    }

    Args args;
};

///////////////////////////////////////////////////////////////////////////////
// operator |
//
template<typename Rng,typename Base,typename Args>
inline BOOST_DEDUCED_TYPENAME Base::BOOST_NESTED_TEMPLATE apply<Rng, Args>::type
operator |(Rng & rng, range_adaptor<Base,Args> const & that)
{
    return that.make_range(rng,that.args);
}

///////////////////////////////////////////////////////////////////////////////
// operator |
//
template<typename Rng,typename Base,typename Args>
inline BOOST_DEDUCED_TYPENAME Base::BOOST_NESTED_TEMPLATE apply<Rng const, Args>::type
operator |(Rng const & rng, range_adaptor<Base,Args> const & that)
{
    return that.make_range(rng,that.args);
}

} // namespace boost

#endif
