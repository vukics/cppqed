///////////////////////////////////////////////////////////////////////////////
/// \file iterator_cast.hpp
///   Defines the iterator_cast() function, as well as the auto_base() helper
//
///////////////////////////////////////////////////////////////////////////////
#ifndef ITERATOR_CAST_EN_12_17_2004_HPP
#define ITERATOR_CAST_EN_12_17_2004_HPP

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>

namespace boost
{

namespace range_ex_detail
{

    /// iterator_cast_result
    template<typename Base,typename Iter>
    struct iterator_cast_result
        : boost::disable_if<boost::is_same<Base,Iter>, Base>
    {
    };

} // namespace range_ex_detail

// work-around for Doxygen bug
#ifndef BOOST_RANGE_EX_DOXYGEN_INVOKED
/// iterator_cast
template<typename Iter>
inline Iter iterator_cast(Iter iter)
{
    return iter;
}

/// iterator_cast
template<typename Base,typename Iter>
inline BOOST_DEDUCED_TYPENAME range_ex_detail::iterator_cast_result<Base, Iter>::type
iterator_cast(Iter iter)
{
    return iterator_cast<Base>(iter.base());
}
#else
/// \brief casts the specified iterator adaptor to the desired base iterator type
///
/// casts the specified iterator adaptor to the desired base iterator type
///
template<typename Base,typename Iter> Base iterator_cast(Iter iter);
#endif

namespace range_ex_detail
{

    /// auto_base_t
    template<typename Iter>
    struct auto_base_t
    {
        explicit auto_base_t(Iter it)
            : iter(it)
        {
        }

        template<typename T>
        operator T() const
        {
            return iterator_cast<T>(iter);
        }

    private:
        Iter iter;
    };

} // namespace range_ex_detail

/// \brief provides an implicit conversion from an iterator adaptor
/// to a base iterator type
///
/// provides an implicit conversion from an iterator adaptor
/// to a base iterator type
///
template<typename Iter>
inline range_ex_detail::auto_base_t<Iter> auto_base(Iter iter)
{
    return range_ex_detail::auto_base_t<Iter>(iter);
}

} // namespace boost

#endif
