///////////////////////////////////////////////////////////////////////////////
/// \file algorithm.hpp
///   Contains range-based versions of the std algorithms
//
/////////////////////////////////////////////////////////////////////////////

// Copyright 2004 Eric Niebler.
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// (C) Copyright Thorsten Ottosen 2002-3. Permission to copy, use, modify,
// sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.

// (C) Copyright Jeremy Siek 2001. Permission to copy, use, modify,
// sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.

// Mutating algorithms originally written by Vladimir Prus'
// <ghost@cs.msu.su> code from Boost Wiki

// Problem: should member functions be called automatically? Or should the user
// know that it is better to call map::find() than find( map )?

#if defined(_MSC_VER) && _MSC_VER >= 1000
#pragma once
#endif

#ifndef ALGORITHM_EN_14_12_2004
#define ALGORITHM_EN_14_12_2004

#include <algorithm>
#include <boost/iterator/iterator_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/range/difference_type.hpp>
#include "./detail/adl_begin_end.hpp"
#include "./detail/has_find.hpp"
#include "./detail/has_remove.hpp"
#include "./detail/has_remove_if.hpp"
#include "./detail/has_unique.hpp"
#include "./detail/has_reverse.hpp"
#include "./detail/has_sort.hpp"
#include "./detail/has_lower_bound.hpp"
#include "./detail/has_upper_bound.hpp"
#include "./detail/has_equal_range.hpp"

namespace boost
{
    namespace range_ex_detail
    {
        template<typename Rng>
        struct iter_pair
        {
            typedef BOOST_DEDUCED_TYPENAME boost::range_result_iterator<
                Rng
            >::type iterator;

            typedef std::pair<iterator,iterator> type;
        };
    }

    /////////////////////////////////////////////////////////////////////////
    // Non-Modifying Sequence Operations
    /////////////////////////////////////////////////////////////////////////

    /// \brief template function for_each
    ///
    /// range-based version of the for_each std algorithm
    ///
    /// \pre Rng meets the requirements for a Input range
    template<typename Rng,typename Fun>
    inline Fun for_each(Rng & rng,Fun fun)
    {
        return std::for_each(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),fun);
    }

    /// \overload
    template<typename Rng,typename Fun>
    inline Fun for_each(Rng const & rng,Fun fun)
    {
        return std::for_each(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),fun);
    }

    namespace range_ex_detail
    {
        template<typename Rng,typename Val>
        inline BOOST_DEDUCED_TYPENAME boost::lazy_enable_if<
            has_find<Rng>
          , boost::range_result_iterator<Rng>
        >::type
        find_impl(Rng & rng,Val const & val)
        {
            return rng.find(val);
        }

        template<typename Rng,typename Val>
        inline BOOST_DEDUCED_TYPENAME boost::lazy_disable_if<
            has_find<Rng>
          , boost::range_result_iterator<Rng>
        >::type
        find_impl(Rng & rng,Val const & val)
        {
            return std::find(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),val);
        }
    }

    /// \brief template function find
    ///
    /// range-based version of the find std algorithm
    ///
    /// \pre Rng meets the requirements for a Input range
    template<typename Rng,typename Val>
    inline BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type
    find(Rng & rng,Val const & val)
    {
        return range_ex_detail::find_impl(rng,val);
    }

    /// \overload
    template<typename Rng,typename Val>
    inline BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type
    find(Rng const & rng,Val const & val)
    {
        return range_ex_detail::find_impl(rng,val);
    }

    /// \brief template function find_if
    ///
    /// range-based version of the find_if std algorithm
    ///
    /// \pre Rng meets the requirements for a Input range
    template<typename Rng,typename Pred>
    inline BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type
    find_if(Rng & rng,Pred pred)
    {
        return std::find_if(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),pred);
    }

    /// \overload
    template<typename Rng,typename Pred>
    inline BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type
    find_if(Rng const & rng,Pred pred)
    {
        return std::find_if(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),pred);
    }

    /// \brief template function find_end
    ///
    /// range-based version of the find_end std algorithm
    ///
    /// \pre Rng1 meets the requirements for a Forward range
    /// \pre Rng2 meets the requirements for a Forward range
    template<typename Rng1,typename Rng2>
    inline BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng1>::type
    find_end(Rng1 & rng1,Rng2 const & rng2)
    {
        return std::find_end(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                             range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2));
    }

    /// \overload
    template<typename Rng1,typename Rng2>
    inline BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng1>::type
    find_end(Rng1 const & rng1,Rng2 const & rng2)
    {
        return std::find_end(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                             range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2));
    }

    /// \overload
    template<typename Rng1,typename Rng2,typename BinPred>
    inline BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng1>::type
    find_end(Rng1 & rng1,Rng2 const & rng2,BinPred pred)
    {
        return std::find_end(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                             range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2),pred);
    }

    /// \overload
    template<typename Rng1,typename Rng2,typename BinPred>
    inline BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng1>::type
    find_end(Rng1 const & rng1,Rng2 const & rng2,BinPred pred)
    {
        return std::find_end(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                             range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2),pred);
    }

    /// \brief template function find_first_of
    ///
    /// range-based version of the find_first_of std algorithm
    ///
    /// \pre Rng1 meets the requirements for a Forward range
    /// \pre Rng2 meets the requirements for a Forward range
    template<typename Rng1,typename Rng2>
    inline BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng1>::type
    find_first_of(Rng1 & rng1,Rng2 const & rng2)
    {
        return std::find_first_of(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                                  range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2));
    }

    /// \overload
    template<typename Rng1,typename Rng2>
    inline BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng1>::type
    find_first_of(Rng1 const & rng1,Rng2 const & rng2)
    {
        return std::find_first_of(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                                  range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2));
    }

    /// \overload
    template<typename Rng1,typename Rng2,typename BinPred>
    inline BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng1>::type
    find_first_of(Rng1 & rng1,Rng2 const & rng2,BinPred pred)
    {
        return std::find_first_of(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                                  range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2),pred);
    }

    /// \overload
    template<typename Rng1,typename Rng2,typename BinPred>
    inline BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng1>::type
    find_first_of(Rng1 const & rng1,Rng2 const & rng2,BinPred pred)
    {
        return std::find_first_of(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                                  range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2),pred);
    }

    /// \brief template function adjacent_find
    ///
    /// range-based version of the adjacent_find std algorithm
    ///
    /// \pre Rng meets the requirements for a Forward range
    template<typename Rng>
    inline BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type
    adjacent_find(Rng & rng)
    {
        return std::adjacent_find(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng>
    inline BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type
    adjacent_find(Rng const & rng)
    {
        return std::adjacent_find(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng,typename BinPred>
    inline BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type
    adjacent_find(Rng & rng,BinPred pred)
    {
        return std::adjacent_find(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),pred);
    }

    /// \overload
    template<typename Rng,typename BinPred>
    inline BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type
    adjacent_find(Rng const & rng,BinPred pred)
    {
        return std::adjacent_find(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),pred);
    }

    /// \brief template function count
    ///
    /// range-based version of the count std algorithm
    ///
    /// \pre Rng meets the requirements for an Input range
    template<typename Rng,typename Val>
    inline BOOST_DEDUCED_TYPENAME boost::range_difference<Rng>::type
    count(Rng & rng,Val const & val)
    {
        return std::count(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),val);
    }

    /// \overload
    template<typename Rng,typename Val>
    inline BOOST_DEDUCED_TYPENAME boost::range_difference<Rng const>::type
    count(Rng const & rng,Val const & val)
    {
        return std::count(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),val);
    }

    /// \brief template function count_if
    ///
    /// range-based version of the count_if std algorithm
    ///
    /// \pre Rng meets the requirements for an Input range
    template<typename Rng,typename Pred>
    inline BOOST_DEDUCED_TYPENAME boost::range_difference<Rng>::type
    count_if(Rng & rng,Pred pred)
    {
        return std::count_if(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),pred);
    }

    /// \overload
    template<typename Rng,typename Pred>
    inline BOOST_DEDUCED_TYPENAME boost::range_difference<Rng const>::type
    count_if(Rng const & rng,Pred pred)
    {
        return std::count_if(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),pred);
    }

    /// \brief template function mismatch
    ///
    /// range-based version of the mismatch std algorithm
    ///
    /// \pre Rng meets the requirements for an Input range
    /// \pre InIter meets the requirements for an Input iterator
    template<typename Rng,typename InIter>
    inline std::pair<BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type,InIter>
    mismatch(Rng & rng,InIter first)
    {
        return std::mismatch(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),first);
    }

    /// \overload
    template<typename Rng,typename InIter>
    inline std::pair<BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type,InIter>
    mismatch(Rng const & rng,InIter first)
    {
        return std::mismatch(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),first);
    }

    /// \overload
    template<typename Rng,typename InIter,typename BinPred>
    inline std::pair<BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type,InIter>
    mismatch(Rng & rng,InIter first,BinPred pred)
    {
        return std::mismatch(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),first,pred);
    }

    /// \overload
    template<typename Rng,typename InIter,typename BinPred>
    inline std::pair<BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type,InIter>
    mismatch(Rng const & rng,InIter first,BinPred pred)
    {
        return std::mismatch(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),first,pred);
    }

    /// \brief template function equal
    ///
    /// range-based version of the equal std algorithm
    ///
    /// \pre Rng meets the requirements for an Input range
    /// \pre InIter meets the requirements for an Input iterator
    template<typename Rng,typename InIter>
    inline bool equal(Rng & rng,InIter first)
    {
        return std::equal(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),first);
    }

    /// \overload
    template<typename Rng,typename InIter>
    inline bool equal(Rng const & rng,InIter first)
    {
        return std::equal(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),first);
    }

    /// \overload
    template<typename Rng,typename InIter,typename BinPred>
    inline bool equal(Rng & rng,InIter first,BinPred pred)
    {
        return std::equal(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),first,pred);
    }

    /// \overload
    template<typename Rng,typename InIter,typename BinPred>
    inline bool equal(Rng const & rng,InIter first,BinPred pred)
    {
        return std::equal(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),first,pred);
    }

    /// \brief template function search
    ///
    /// range-based version of the search std algorithm
    ///
    /// \pre Rng1 meets the requirements for a Forward range
    /// \pre Rng2 meets the requirements for a Forward range
    template<typename Rng1,typename Rng2>
    inline BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng1>::type
    search(Rng1 & rng1,Rng2 const & rng2)
    {
        return std::search(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                           range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2));
    }

    /// \overload
    template<typename Rng1,typename Rng2>
    inline BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng1>::type
    search(Rng1 const & rng1,Rng2 const & rng2)
    {
        return std::search(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                           range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2));
    }

    /// \overload
    template<typename Rng1,typename Rng2,typename BinPred>
    inline BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng1>::type
    search(Rng1 & rng1,Rng2 const & rng2,BinPred pred)
    {
        return std::search(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                           range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2),pred);
    }

    /// \overload
    template<typename Rng1,typename Rng2,typename BinPred>
    inline BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng1>::type
    search(Rng1 const & rng1,Rng2 const & rng2,BinPred pred)
    {
        return std::search(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                           range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2),pred);
    }

    /////////////////////////////////////////////////////////////////////////
    // Modifying Sequence Operations
    /////////////////////////////////////////////////////////////////////////

    /// \brief template function copy
    ///
    /// range-based version of the copy std algorithm
    ///
    /// \pre Rng meets the requirements for an Input range
    /// \pre OutIter meets the requirements for an Output iterator
    template<typename Rng,typename OutIter>
    inline OutIter copy(Rng const & rng,OutIter out)
    {
        return std::copy(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),out);
    }

    /// \brief template function copy_backwards
    ///
    /// range-based version of the copy_backwards std algorithm
    ///
    /// \pre Rng meets the requirements for a Bidirectional range
    /// \pre BidiIter meets the requirements for a Bidirectional iterator
    template<typename Rng,typename BidiIter>
    inline BidiIter copy_backward(Rng const & rng,BidiIter out)
    {
        return std::copy_backward(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),out);
    }

    /// \brief template function transform
    ///
    /// range-based version of the transform std algorithm
    ///
    /// \pre Rng meets the requirements for an Input range
    /// \pre InIter meets the requirements for an Input iterator
    /// \pre OutIter meets the requirements for an Output iterator
    template<typename Rng,typename OutIter,typename UnaryOp>
    inline OutIter transform(Rng const & rng,OutIter out,UnaryOp fun)
    {
        return std::transform(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),out,fun);
    }

    /// \overload
    template<typename Rng,typename InIter,typename OutIter,typename BinOp>
    inline OutIter transform(Rng const & rng,InIter first2,OutIter out,BinOp fun)
    {
        return std::transform(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),first2,out,fun);
    }

    /// \brief template function replace
    ///
    /// range-based version of the replace std algorithm
    ///
    /// \pre Rng meets the requirements for a Forward range
    template<typename Rng,typename Val>
    inline void replace(Rng & rng,Val const & what,Val const & with_what)
    {
        return std::replace(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),what,with_what);
    }

    /// \overload
    template<typename Rng,typename Val>
    inline void replace(Rng const & rng,Val const & what,Val const & with_what)
    {
        return std::replace(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),what,with_what);
    }

    /// \brief template function replace_if
    ///
    /// range-based version of the replace_if std algorithm
    ///
    /// \pre Rng meets the requirements for a Forward range
    template<typename Rng,typename Pred,typename Val>
    inline void replace_if(Rng & rng,Pred pred,Val const & val)
    {
        return std::replace_if(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),pred,val);
    }

    /// \overload
    template<typename Rng,typename Pred,typename Val>
    inline void replace_if(Rng const & rng,Pred pred,Val const & val)
    {
        return std::replace_if(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),pred,val);
    }

    /// \brief template function replace_copy
    ///
    /// range-based version of the replace_copy std algorithm
    ///
    /// \pre Rng meets the requirements for an Input range
    /// \pre OutIter meets the requirements for an Output iterator
    template<typename Rng,typename OutIter,typename Val>
    inline OutIter replace_copy(Rng const & rng,OutIter out,Val const & what,Val const & with_what)
    {
        return std::replace_copy(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),out,what,with_what);
    }

    /// \brief template function replace_copy_if
    ///
    /// range-based version of the replace_copy_if std algorithm
    ///
    /// \pre Rng meets the requirements for an Input range
    /// \pre OutIter meets the requirements for an Output iterator
    template<typename Rng,typename OutIter,typename Pred,typename Val>
    inline OutIter replace_copy_if(Rng const & rng,OutIter out,Pred pred,Val const & val)
    {
        return std::replace_copy_if(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),out,pred,val);
    }

    /// \brief template function fill
    ///
    /// range-based version of the fill std algorithm
    ///
    /// \pre Rng meets the requirements for a Forward range
    template<typename Rng,typename Val>
    inline void fill(Rng & rng,Val const & val)
    {
        std::fill(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),val);
    }

    /// \overload
    template<typename Rng,typename Val>
    inline void fill(Rng const & rng,Val const & val)
    {
        std::fill(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),val);
    }

    /// \brief template function fill_n
    ///
    /// range-based version of the fill_n std algorithm
    ///
    /// \pre Rng meets the requirements for an Output range
    template<typename Rng,typename Int,typename Val>
    inline void fill_n(Rng & rng,Int size,Val const & val)
    {
        // BUGBUG an Output range? Rethink this, and output ranges in general.
        std::fill_n(range_ex_detail::adl_begin(rng),size,val);
    }

    /// \overload
    template<typename Rng,typename Int,typename Val>
    inline void fill_n(Rng const & rng,Int size,Val const & val)
    {
        std::fill_n(range_ex_detail::adl_begin(rng),size,val);
    }

    /// \brief template function generate
    ///
    /// range-based version of the generate std algorithm
    ///
    /// \pre Rng meets the requirements for a Forward range
    template<typename Rng,typename Generator>
    inline void generate(Rng & rng,Generator gen)
    {
        std::generate(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),gen);
    }

    /// \overload
    template<typename Rng,typename Generator>
    inline void generate(Rng const & rng,Generator gen)
    {
        std::generate(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),gen);
    }

    /// \brief template function generate_n
    ///
    /// range-based version of the generate_n std algorithm
    ///
    /// \pre Rng meets the requirements for an Output range
    template<typename Rng,typename Int,typename Generator>
    void generate_n(Rng & rng,Int size,Generator gen)
    {
        std::generate_n(range_ex_detail::adl_begin(rng),size,gen);
    }

    /// \overload
    template<typename Rng,typename Int,typename Generator>
    void generate_n(Rng const & rng,Int size,Generator gen)
    {
        std::generate_n(range_ex_detail::adl_begin(rng),size,gen);
    }

    namespace range_ex_detail
    {
        template<typename Rng,typename Val>
        inline BOOST_DEDUCED_TYPENAME boost::lazy_enable_if<
            has_remove<Rng>
          , boost::range_result_iterator<Rng>
        >::type
        remove_impl(Rng & rng,Val const & val)
        {
            rng.remove(val);
            return range_ex_detail::adl_end(rng);
        }

        template<typename Rng,typename Val>
        inline BOOST_DEDUCED_TYPENAME boost::lazy_disable_if<
            has_remove<Rng>
          , boost::range_result_iterator<Rng>
        >::type
        remove_impl(Rng & rng,Val const & val)
        {
            return std::remove(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),val);
        }
    }

    /// \brief template function remove
    ///
    /// range-based version of the remove std algorithm
    ///
    /// \pre Rng meets the requirements for a Forward range
    template<typename Rng,typename Val>
    inline BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type
    remove(Rng & rng,Val const & val)
    {
        return range_ex_detail::remove_impl(rng,val);
    }

    /// \overload
    template<typename Rng,typename Val>
    inline BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type
    remove(Rng const & rng,Val const & val)
    {
        return range_ex_detail::remove_impl(rng,val);
    }

    namespace range_ex_detail
    {
        template<typename Rng,typename Pred>
        inline BOOST_DEDUCED_TYPENAME boost::lazy_enable_if<
            has_remove_if<Rng>
          , boost::range_result_iterator<Rng>
        >::type
        remove_if_impl(Rng & rng,Pred pred)
        {
            rng.remove_if(pred);
            return range_ex_detail::adl_end(rng);
        }

        template<typename Rng,typename Pred>
        inline BOOST_DEDUCED_TYPENAME boost::lazy_disable_if<
            has_remove_if<Rng>
          , boost::range_result_iterator<Rng>
        >::type
        remove_if_impl(Rng & rng,Pred pred)
        {
            return std::remove_if(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),pred);
        }
    }

    /// \brief template function remove_if
    ///
    /// range-based version of the remove_if std algorithm
    ///
    /// \pre Rng meets the requirements for a Forward range
    template<typename Rng,typename Pred>
    inline BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type
    remove_if(Rng & rng,Pred pred)
    {
        return range_ex_detail::remove_if_impl(rng,pred);
    }

    /// \overload
    template<typename Rng,typename Pred>
    inline BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type
    remove_if(Rng const & rng,Pred pred)
    {
        return range_ex_detail::remove_if_impl(rng,pred);
    }

    /// \brief template function remove_copy
    ///
    /// range-based version of the remove_copy std algorithm
    ///
    /// \pre Rng meets the requirements for an Input range
    /// \pre OutIter meets the requirements for an Output iterator
    template<typename Rng,typename OutIter,typename Val>
    inline OutIter remove_copy(Rng const & rng,OutIter out,Val const & val)
    {
        return std::remove_copy(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),out,val);
    }

    /// \brief template function remove_copy_if
    ///
    /// range-based version of the remove_copy_if std algorithm
    ///
    /// \pre Rng meets the requirements for an Input range
    /// \pre OutIter meets the requirements for an Output iterator
    template<typename Rng,typename OutIter,typename Pred>
    inline OutIter remove_copy_if(Rng const & rng,OutIter out,Pred pred)
    {
        return std::remove_copy_if(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),out,pred);
    }

    namespace range_ex_detail
    {
        template<typename Rng>
        inline BOOST_DEDUCED_TYPENAME boost::lazy_enable_if<
            has_unique<Rng>
          , boost::range_result_iterator<Rng>
        >::type
        unique_impl(Rng & rng)
        {
            rng.unique();
            return range_ex_detail::adl_end(rng);
        }

        template<typename Rng>
        inline BOOST_DEDUCED_TYPENAME boost::lazy_disable_if<
            has_unique<Rng>
          , boost::range_result_iterator<Rng>
        >::type
        unique_impl(Rng & rng)
        {
            return std::unique(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng));
        }

        template<typename Rng,typename Pred>
        inline BOOST_DEDUCED_TYPENAME boost::lazy_enable_if<
            has_unique<Rng>
          , boost::range_result_iterator<Rng>
        >::type
        unique_if_impl(Rng & rng,Pred pred)
        {
            rng.unique(pred);
            return range_ex_detail::adl_end(rng);
        }

        template<typename Rng,typename Pred>
        inline BOOST_DEDUCED_TYPENAME boost::lazy_disable_if<
            has_unique<Rng>
          , boost::range_result_iterator<Rng>
        >::type
        unique_if_impl(Rng & rng,Pred pred)
        {
            return std::unique(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),pred);
        }
    }

    /// \brief template function unique
    ///
    /// range-based version of the unique std algorithm
    ///
    /// \pre Rng meets the requirements for a Forward range
    template<typename Rng>
    inline BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type
    unique(Rng & rng)
    {
        return range_ex_detail::unique_impl(rng);
    }

    /// \overload
    template<typename Rng>
    inline BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type
    unique(Rng const & rng)
    {
        return range_ex_detail::unique_impl(rng);
    }

    /// \overload
    template<typename Rng,typename Pred>
    inline BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type
    unique(Rng & rng,Pred pred)
    {
        return range_ex_detail::unique_if_impl(rng,pred);
    }

    /// \overload
    template<typename Rng,typename Pred>
    inline BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type
    unique(Rng const & rng,Pred pred)
    {
        return range_ex_detail::unique_if_impl(rng,pred);
    }

    /// \brief template function unique_copy
    ///
    /// range-based version of the unique_copy std algorithm
    ///
    /// \pre Rng meets the requirements for an Input range
    /// \pre OutIter meets the requirements for an Output iterator
    template<typename Rng,typename OutIter>
    inline OutIter unique_copy(Rng const & rng,OutIter out)
    {
        return std::unique_copy(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),out);
    }

    /// \overload
    template<typename Rng,typename OutIter,typename Pred>
    inline OutIter unique_copy(Rng const & rng,OutIter out,Pred pred)
    {
        return std::unique_copy(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),out,pred);
    }

    namespace range_ex_detail
    {
        template<typename Rng>
        inline BOOST_DEDUCED_TYPENAME boost::enable_if<
            has_reverse<Rng>
        >::type
        reverse_impl(Rng & rng)
        {
            rng.reverse();
        }

        template<typename Rng>
        inline BOOST_DEDUCED_TYPENAME boost::disable_if<
            has_reverse<Rng>
        >::type
        reverse_impl(Rng & rng)
        {
            std::reverse(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng));
        }
    }

    /// \brief template function reverse
    ///
    /// range-based version of the reverse std algorithm
    ///
    /// \pre Rng meets the requirements for a Bidirectional range
    template<typename Rng>
    inline void reverse(Rng & rng)
    {
        range_ex_detail::reverse_impl(rng);
    }

    /// \overload
    template<typename Rng>
    inline void reverse(Rng const & rng)
    {
        range_ex_detail::reverse_impl(rng);
    }

    /// \brief template function reverse_copy
    ///
    /// range-based version of the reverse_copy std algorithm
    ///
    /// \pre Rng meets the requirements for a Bidirectional range
    /// \pre OutIter meets the requirements for an Output iterator
    template<typename Rng,typename OutIter>
    inline OutIter reverse_copy(Rng const & rng,OutIter out)
    {
        return std::reverse_copy(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),out);
    }

    /// \brief template function rotate
    ///
    /// range-based version of the rotate std algorithm
    ///
    /// \pre Rng meets the requirements for a Forward range
    template<typename Rng>
    inline void rotate(Rng & rng,
        BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type middle)
    {
        std::rotate(range_ex_detail::adl_begin(rng),middle,range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng>
    inline void rotate(Rng const & rng,
        BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type middle)
    {
        std::rotate(range_ex_detail::adl_begin(rng),middle,range_ex_detail::adl_end(rng));
    }

    /// \brief template function rotate_copy
    ///
    /// range-based version of the rotate_copy std algorithm
    ///
    /// \pre Rng meets the requirements for a Forward range
    /// \pre OutIter meets the requirements for an Output iterator
    template<typename Rng,typename OutIter>
    inline OutIter rotate_copy(Rng & rng,
        BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type middle,
        OutIter out)
    {
        return std::rotate_copy(range_ex_detail::adl_begin(rng),middle,range_ex_detail::adl_end(rng),out);
    }

    /// \overload
    template<typename Rng,typename OutIter>
    inline OutIter rotate_copy(Rng const & rng,
        BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type middle,
        OutIter out)
    {
        return std::rotate_copy(range_ex_detail::adl_begin(rng),middle,range_ex_detail::adl_end(rng),out);
    }

    /// \brief template function random_shuffle
    ///
    /// range-based version of the random_shuffle std algorithm
    ///
    /// \pre Rng meets the requirements for a Random Access range
    template<typename Rng>
    inline void random_shuffle(Rng & rng)
    {
        std::random_shuffle(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng>
    inline void random_shuffle(Rng const & rng)
    {
        std::random_shuffle(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng,typename Generator>
    inline void random_shuffle(Rng & rng,Generator gen)
    {
        std::random_shuffle(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),gen);
    }

    /// \overload
    template<typename Rng,typename Generator>
    inline void random_shuffle(Rng const & rng,Generator gen)
    {
        std::random_shuffle(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),gen);
    }

    /// \brief template function partition
    ///
    /// range-based version of the partition std algorithm
    ///
    /// \pre Rng meets the requirements for a Bidirectional range
    template<typename Rng,typename Pred>
    inline BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type
    partition(Rng & rng,Pred pred)
    {
        return std::partition(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),pred);
    }

    /// \overload
    template<typename Rng,typename Pred>
    inline BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type
    partition(Rng const & rng,Pred pred)
    {
        return std::partition(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),pred);
    }

    /// \brief template function stable_partition
    ///
    /// range-based version of the stable_partition std algorithm
    ///
    /// \pre Rng meets the requirements for a Bidirectional range
    template<typename Rng,typename Pred>
    inline BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type
    stable_partition(Rng & rng,Pred pred)
    {
        return std::stable_partition(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),pred);
    }

    /// \overload
    template<typename Rng,typename Pred>
    inline BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type
    stable_partition(Rng const & rng,Pred pred)
    {
        return std::stable_partition(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),pred);
    }

    namespace range_ex_detail
    {
        template<typename Rng>
        inline BOOST_DEDUCED_TYPENAME boost::enable_if<
            has_sort<Rng>
        >::type
        sort_impl(Rng & rng)
        {
            rng.sort();
        }

        template<typename Rng>
        inline BOOST_DEDUCED_TYPENAME boost::disable_if<
            has_sort<Rng>
        >::type
        sort_impl(Rng & rng)
        {
            std::sort(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng));
        }
    }

    /// \brief template function sort
    ///
    /// range-based version of the sort std algorithm
    ///
    /// \pre Rng meets the requirements for a Random Access range
    template<typename Rng>
    inline void sort(Rng & rng)
    {
        range_ex_detail::sort_impl(rng);
    }

    /// \overload
    template<typename Rng>
    inline void sort(Rng const & rng)
    {
        range_ex_detail::sort_impl(rng);
    }

    namespace range_ex_detail
    {
        template<typename Rng,typename Cmp>
        inline BOOST_DEDUCED_TYPENAME boost::enable_if<
            has_sort<Rng>
        >::type
        sort_impl(Rng & rng,Cmp cmp)
        {
            rng.sort(cmp);
        }

        template<typename Rng,typename Cmp>
        inline BOOST_DEDUCED_TYPENAME boost::disable_if<
            has_sort<Rng>
        >::type
        sort_impl(Rng & rng,Cmp cmp)
        {
            std::sort(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),cmp);
        }
    }

    /// \overload
    template<typename Rng,typename Cmp>
    inline void sort(Rng & rng,Cmp cmp)
    {
        range_ex_detail::sort_impl(rng,cmp);
    }

    /// \overload
    template<typename Rng,typename Cmp>
    inline void sort(Rng const & rng,Cmp cmp)
    {
        range_ex_detail::sort_impl(rng,cmp);
    }

    /// \brief template function stable_sort
    ///
    /// range-based version of the stable_sort std algorithm
    ///
    /// \pre Rng meets the requirements for a Random Access range
    template<typename Rng>
    inline void stable_sort(Rng & rng)
    {
        std::stable_sort(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng>
    inline void stable_sort(Rng const & rng)
    {
        std::stable_sort(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng,typename Cmp>
    inline void stable_sort(Rng & rng,Cmp cmp)
    {
        std::stable_sort(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),cmp);
    }

    /// \overload
    template<typename Rng,typename Cmp>
    inline void stable_sort(Rng const & rng,Cmp cmp)
    {
        std::stable_sort(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),cmp);
    }

    /// \brief template function partial_sort
    ///
    /// range-based version of the partial_sort std algorithm
    ///
    /// \pre Rng meets the requirements for a Random Access range
    template<typename Rng>
    inline void partial_sort(Rng & rng,
        BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type middle)
    {
        std::partial_sort(range_ex_detail::adl_begin(rng),middle,range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng>
    inline void partial_sort(Rng const & rng,
        BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type middle)
    {
        std::partial_sort(range_ex_detail::adl_begin(rng),middle,range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng,typename Cmp>
    inline void partial_sort(Rng & rng,
        BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type middle,
        Cmp cmp)
    {
        std::partial_sort(range_ex_detail::adl_begin(rng),middle,range_ex_detail::adl_end(rng),cmp);
    }

    /// \overload
    template<typename Rng,typename Cmp>
    inline void partial_sort(Rng const & rng,
        BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type middle,
        Cmp cmp)
    {
        std::partial_sort(range_ex_detail::adl_begin(rng),middle,range_ex_detail::adl_end(rng),cmp);
    }

    /// \brief template function partial_sort_copy
    ///
    /// range-based version of the partial_sort_copy std algorithm
    ///
    /// \pre Rng1 meets the requirements for a Input range
    /// \pre Rng2 meets the requirements for a Random Access range
    template<typename Rng1,typename Rng2>
    inline BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng2>::type
    partial_sort_copy(Rng1 const & rng1,Rng2 & rng2)
    {
        return std::partial_sort_copy(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                                      range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2));
    }

    /// \overload
    template<typename Rng1,typename Rng2>
    inline BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng2>::type
    partial_sort_copy(Rng1 const & rng1,Rng2 const & rng2)
    {
        return std::partial_sort_copy(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                                      range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2));
    }

    /// \overload
    template<typename Rng1,typename Rng2,typename Cmp>
    inline BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng2>::type
    partial_sort_copy(Rng1 const & rng1,Rng2 & rng2,Cmp cmp)
    {
        return std::partial_sort_copy(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                                      range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2),cmp);
    }

    /// \overload
    template<typename Rng1,typename Rng2,typename Cmp>
    inline BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng2>::type
    partial_sort_copy(Rng1 const & rng1,Rng2 const & rng2,Cmp cmp)
    {
        return std::partial_sort_copy(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                                      range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2),cmp);
    }

    /// \brief template function nth_element
    ///
    /// range-based version of the nth_element std algorithm
    ///
    /// \pre Rng meets the requirements for a Random Access range
    template<typename Rng>
    inline void nth_element(Rng & rng,
        BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type nth)
    {
        std::nth_element(range_ex_detail::adl_begin(rng),nth,range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng>
    inline void nth_element(Rng const & rng,
        BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type nth)
    {
        std::nth_element(range_ex_detail::adl_begin(rng),nth,range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng,typename Cmp>
    inline void nth_element(Rng & rng,
        BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type nth,
        Cmp cmp)
    {
        std::nth_element(range_ex_detail::adl_begin(rng),nth,range_ex_detail::adl_end(rng),cmp);
    }

    /// \overload
    template<typename Rng,typename Cmp>
    inline void nth_element(Rng const & rng,
        BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type nth,
        Cmp cmp)
    {
        std::nth_element(range_ex_detail::adl_begin(rng),nth,range_ex_detail::adl_end(rng),cmp);
    }

    namespace range_ex_detail
    {
        template<typename Rng,typename Val>
        inline BOOST_DEDUCED_TYPENAME boost::lazy_enable_if<
            has_lower_bound<Rng>
          , boost::range_result_iterator<Rng>
        >::type
        lower_bound_impl(Rng & rng,Val const & val)
        {
            return rng.lower_bound(val);
        }

        template<typename Rng,typename Val>
        inline BOOST_DEDUCED_TYPENAME boost::lazy_disable_if<
            has_lower_bound<Rng>
          , boost::range_result_iterator<Rng>
        >::type
        lower_bound_impl(Rng & rng,Val const & val)
        {
            return std::lower_bound(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),val);
        }
    }

    /// \brief template function lower_bound
    ///
    /// range-based version of the lower_bound std algorithm
    ///
    /// \pre Rng meets the requirements for a Forward range
    template<typename Rng,typename Val>
    inline BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type
    lower_bound(Rng & rng,Val const & val)
    {
        return range_ex_detail::lower_bound_impl(rng,val);
    }

    /// \overload
    template<typename Rng,typename Val>
    inline BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type
    lower_bound(Rng const & rng,Val const & val)
    {
        return range_ex_detail::lower_bound_impl(rng,val);
    }

    /// \overload
    template<typename Rng,typename Val,typename Cmp>
    inline BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type
    lower_bound(Rng & rng,Val const & val,Cmp cmp)
    {
        return std::lower_bound(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),val,cmp);
    }

    /// \overload
    template<typename Rng,typename Val,typename Cmp>
    inline BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type
    lower_bound(Rng const & rng,Val const & val,Cmp cmp)
    {
        return std::lower_bound(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),val,cmp);
    }

    namespace range_ex_detail
    {
        template<typename Rng,typename Val>
        inline BOOST_DEDUCED_TYPENAME boost::lazy_enable_if<
            has_upper_bound<Rng>
          , boost::range_result_iterator<Rng>
        >::type
        upper_bound_impl(Rng & rng,Val const & val)
        {
            return rng.upper_bound(val);
        }

        template<typename Rng,typename Val>
        inline BOOST_DEDUCED_TYPENAME boost::lazy_disable_if<
            has_upper_bound<Rng>
          , boost::range_result_iterator<Rng>
        >::type
        upper_bound_impl(Rng & rng,Val const & val)
        {
            return std::upper_bound(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),val);
        }
    }

    /// \brief template function upper_bound
    ///
    /// range-based version of the upper_bound std algorithm
    ///
    /// \pre Rng meets the requirements for a Forward range
    template<typename Rng,typename Val>
    inline BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type
    upper_bound(Rng & rng,Val const & val)
    {
        return range_ex_detail::upper_bound_impl(rng,val);
    }

    /// \overload
    template<typename Rng,typename Val>
    inline BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type
    upper_bound(Rng const & rng,Val const & val)
    {
        return range_ex_detail::upper_bound_impl(rng,val);
    }

    /// \overload
    template<typename Rng,typename Val,typename Cmp>
    inline BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type
    upper_bound(Rng & rng,Val const & val,Cmp cmp)
    {
        return std::upper_bound(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),val,cmp);
    }

    /// \overload
    template<typename Rng,typename Val,typename Cmp>
    inline BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type
    upper_bound(Rng const & rng,Val const & val,Cmp cmp)
    {
        return std::upper_bound(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),val,cmp);
    }

    namespace range_ex_detail
    {
        template<typename Rng,typename Val>
        inline BOOST_DEDUCED_TYPENAME boost::lazy_enable_if<
            has_equal_range<Rng>
          , iter_pair<Rng>
        >::type
        equal_range_impl(Rng & rng,Val const & val)
        {
            return rng.equal_range(val);
        }

        template<typename Rng,typename Val>
        inline BOOST_DEDUCED_TYPENAME boost::lazy_disable_if<
            has_equal_range<Rng>
          , iter_pair<Rng>
        >::type
        equal_range_impl(Rng & rng,Val const & val)
        {
            return std::equal_range(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),val);
        }
    }

    /// \brief template function equal_range
    ///
    /// range-based version of the equal_range std algorithm
    ///
    /// \pre Rng meets the requirements for a Forward range
    template<typename Rng,typename Val>
    inline BOOST_DEDUCED_TYPENAME range_ex_detail::iter_pair<Rng>::type
    equal_range(Rng & rng,Val const & val)
    {
        return range_ex_detail::equal_range_impl(rng,val);
    }

    /// \overload
    template<typename Rng,typename Val>
    inline BOOST_DEDUCED_TYPENAME range_ex_detail::iter_pair<Rng const>::type
    equal_range(Rng const & rng,Val const & val)
    {
        return range_ex_detail::equal_range_impl(rng,val);
    }

    /// \overload
    template<typename Rng,typename Val,typename Cmp>
    inline BOOST_DEDUCED_TYPENAME range_ex_detail::iter_pair<Rng>::type
    equal_range(Rng & rng,Val const & val,Cmp cmp)
    {
        return std::equal_range(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),val,cmp);
    }

    /// \overload
    template<typename Rng,typename Val,typename Cmp>
    inline BOOST_DEDUCED_TYPENAME range_ex_detail::iter_pair<Rng const>::type
    equal_range(Rng const & rng,Val const & val,Cmp cmp)
    {
        return std::equal_range(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),val,cmp);
    }

    /// \brief template function binary_search
    ///
    /// range-based version of the binary_search std algorithm
    ///
    /// \pre Rng meets the requirements for a Forward range
    template<typename Rng,typename Val>
    inline bool binary_search(Rng const & rng,Val const & val)
    {
        return std::binary_search(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),val);
    }

    /// \overload
    template<typename Rng,typename Val,typename Cmp>
    inline bool binary_search(Rng const & rng,Val const & val,Cmp cmp)
    {
        return std::binary_search(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),val,cmp);
    }

    /// \brief template function merge
    ///
    /// range-based version of the merge std algorithm
    ///
    template<typename Rng1,typename Rng2,typename OutIter>
    inline OutIter merge(Rng1 const & rng1,Rng2 const & rng2,OutIter out)
    {
        return std::merge(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                          range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2),out);
    }

    /// \overload
    template<typename Rng1,typename Rng2,typename OutIter,typename Cmp>
    inline OutIter merge(Rng1 const & rng1,Rng2 const & rng2,OutIter out,Cmp cmp)
    {
        return std::merge(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                          range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2),out,cmp);
    }

    /// \brief template function inplace_merge
    ///
    /// range-based version of the inplace_merge std algorithm
    ///
    /// \pre Rng meets the requirements for a Bidirectional range
    template<typename Rng>
    inline void inplace_merge(Rng & rng,
        BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type middle)
    {
        std::inplace_merge(range_ex_detail::adl_begin(rng),middle,range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng>
    inline void inplace_merge(Rng const & rng,
        BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type middle)
    {
        std::inplace_merge(range_ex_detail::adl_begin(rng),middle,range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng,typename Cmp>
    inline void inplace_merge(Rng & rng,
        BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type middle,
        Cmp cmp)
    {
        std::inplace_merge(range_ex_detail::adl_begin(rng),middle,range_ex_detail::adl_end(rng),cmp);
    }

    /// \overload
    template<typename Rng,typename Cmp>
    inline void inplace_merge(Rng const & rng,
        BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type middle,
        Cmp cmp)
    {
        std::inplace_merge(range_ex_detail::adl_begin(rng),middle,range_ex_detail::adl_end(rng),cmp);
    }

    /////////////////////////////////////////////////////////////////////////
    // Set Algorithms
    /////////////////////////////////////////////////////////////////////////

    /// \brief template function includes
    ///
    /// range-based version of the includes std algorithm
    ///
    /// \pre Rng1 meets the requirements for an Input range
    /// \pre Rng2 meets the requirements for an Input range
    template<typename Rng1,typename Rng2>
    inline bool includes(Rng1 const & rng1,Rng2 const & rng2)
    {
        return std::includes(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                             range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2));
    }

    /// \overload
    template<typename Rng1,typename Rng2,typename Cmp>
    inline bool includes(Rng1 const & rng1,Rng2 const & rng2,Cmp cmp)
    {
        return std::includes(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                            range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2),cmp);
    }

    /// \brief template function set_union
    ///
    /// range-based version of the set_union std algorithm
    ///
    /// \pre Rng1 meets the requirements for an Input range
    /// \pre Rng2 meets the requirements for an Input range
    template<typename Rng1,typename Rng2,typename OutIter>
    inline OutIter set_union(Rng1 const & rng1,Rng2 const & rng2,OutIter out)
    {
        return std::set_union(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                              range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2),out);
    }

    /// \overload
    template<typename Rng1,typename Rng2,typename OutIter,typename Cmp>
    inline OutIter set_union(Rng1 const & rng1,Rng2 const & rng2,OutIter out,Cmp cmp)
    {
        return std::set_union(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                              range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2),out,cmp);
    }

    /// \brief template function set_intersection
    ///
    /// range-based version of the set_intersection std algorithm
    ///
    /// \pre Rng1 meets the requirements for an Input range
    /// \pre Rng2 meets the requirements for an Input range
    template<typename Rng1,typename Rng2,typename OutIter>
    inline OutIter set_intersection(Rng1 const & rng1,Rng2 const & rng2,OutIter out)
    {
        return std::set_intersection(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                                     range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2),out);
    }

    /// \overload
    template<typename Rng1,typename Rng2,typename OutIter,typename Cmp>
    inline OutIter set_intersection(Rng1 const & rng1,Rng2 const & rng2,OutIter out,Cmp cmp)
    {
        return std::set_intersection(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                                     range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2),out,cmp);
    }

    /// \brief template function set_difference
    ///
    /// range-based version of the set_difference std algorithm
    ///
    /// \pre Rng1 meets the requirements for an Input range
    /// \pre Rng2 meets the requirements for an Input range
    template<typename Rng1,typename Rng2,typename OutIter>
    inline OutIter set_difference(Rng1 const & rng1,Rng2 const & rng2,OutIter out)
    {
        return std::set_difference(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                                   range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2),out);
    }

    /// \overload
    template<typename Rng1,typename Rng2,typename OutIter,typename Cmp>
    inline OutIter set_difference(Rng1 const & rng1,Rng2 const & rng2,OutIter out,Cmp cmp)
    {
        return std::set_difference(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                                   range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2),out,cmp);
    }

    /// \brief template function set_symmetric_difference
    ///
    /// range-based version of the set_symmetric_difference std algorithm
    ///
    /// \pre Rng1 meets the requirements for an Input range
    /// \pre Rng2 meets the requirements for an Input range
    template<typename Rng1,typename Rng2,typename OutIter>
    inline OutIter set_symmetric_difference(Rng1 const & rng1,Rng2 const & rng2,OutIter out)
    {
        return std::set_symmetric_difference(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                                             range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2),out);
    }

    /// \overload
    template<typename Rng1,typename Rng2,typename OutIter,typename Cmp>
    inline OutIter set_symmetric_difference(Rng1 const & rng1,Rng2 const & rng2,OutIter out,Cmp cmp)
    {
        return std::set_symmetric_difference(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                                             range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2),out,cmp);
    }

    ///////////////////////////////////////////////////////////////////////////
    // Heap Operations
    ///////////////////////////////////////////////////////////////////////////

    /// \brief template function push_heap
    ///
    /// range-based version of the push_heap std algorithm
    ///
    /// \pre Rng meets the requirements for a Random Access range
    template<typename Rng>
    inline void push_heap(Rng & rng)
    {
        std::push_heap(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng>
    inline void push_heap(Rng const & rng)
    {
        std::push_heap(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng,typename Cmp>
    inline void push_heap(Rng & rng,Cmp cmp)
    {
        std::push_heap(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),cmp);
    }

    /// \overload
    template<typename Rng,typename Cmp>
    inline void push_heap(Rng const & rng,Cmp cmp)
    {
        std::push_heap(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),cmp);
    }

    /// \brief template function pop_heap
    ///
    /// range-based version of the pop_heap std algorithm
    ///
    /// \pre Rng meets the requirements for a Random Access range
    template<typename Rng>
    inline void pop_heap(Rng & rng)
    {
        std::pop_heap(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng>
    inline void pop_heap(Rng const & rng)
    {
        std::pop_heap(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng,typename Cmp>
    inline void pop_heap(Rng & rng,Cmp cmp)
    {
        std::pop_heap(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),cmp);
    }

    /// \overload
    template<typename Rng,typename Cmp>
    inline void pop_heap(Rng const & rng,Cmp cmp)
    {
        std::pop_heap(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),cmp);
    }

    /// \brief template function make_heap
    ///
    /// range-based version of the make_heap std algorithm
    ///
    /// \pre Rng meets the requirements for a Random Access range
    template<typename Rng>
    inline void make_heap(Rng & rng)
    {
        std::make_heap(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng>
    inline void make_heap(Rng const & rng)
    {
        std::make_heap(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng,typename Cmp>
    inline void make_heap(Rng & rng,Cmp cmp)
    {
        std::make_heap(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),cmp);
    }

    /// \overload
    template<typename Rng,typename Cmp>
    inline void make_heap(Rng const & rng,Cmp cmp)
    {
        std::make_heap(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),cmp);
    }

    /// \brief template function sort_heap
    ///
    /// range-based version of the sort_heap std algorithm
    ///
    template<typename Rng>
    inline void sort_heap(Rng & rng)
    {
        std::sort_heap(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng>
    inline void sort_heap(Rng const & rng)
    {
        std::sort_heap(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng,typename Cmp>
    inline void sort_heap(Rng & rng,Cmp cmp)
    {
        std::sort_heap(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),cmp);
    }

    /// \overload
    template<typename Rng,typename Cmp>
    inline void sort_heap(Rng const & rng,Cmp cmp)
    {
        std::sort_heap(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),cmp);
    }

    /////////////////////////////////////////////////////////////////////////
    // Minimum and Maximum
    /////////////////////////////////////////////////////////////////////////

    /// \brief template function min_element
    ///
    /// range-based version of the min_element std algorithm
    ///
    /// \pre Rng meets the requirements for a Forward range
    template<typename Rng>
    inline BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type
    min_element(Rng & rng)
    {
        return std::min_element(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng>
    inline BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type
    min_element(Rng const & rng)
    {
        return std::min_element(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng,typename BinPred>
    inline BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type
    min_element(Rng & rng,BinPred pred)
    {
        return std::min_element(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),pred);
    }

    /// \overload
    template<typename Rng,typename BinPred>
    inline BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type
    min_element(Rng const & rng,BinPred pred)
    {
        return std::min_element(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),pred);
    }

    /// \brief template function max_element
    ///
    /// range-based version of the max_element std algorithm
    ///
    /// \pre Rng meets the requirements for a Forward range
    template<typename Rng>
    inline BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type
    max_element(Rng & rng)
    {
        return std::max_element(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng>
    inline BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type
    max_element(Rng const & rng)
    {
        return std::max_element(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng,typename BinPred>
    inline BOOST_DEDUCED_TYPENAME boost::range_iterator<Rng>::type
    max_element(Rng & rng,BinPred pred)
    {
        return std::max_element(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),pred);
    }

    /// \overload
    template<typename Rng,typename BinPred>
    inline BOOST_DEDUCED_TYPENAME boost::range_const_iterator<Rng>::type
    max_element(Rng const & rng,BinPred pred)
    {
        return std::max_element(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),pred);
    }

    /// \brief template function lexicographic_compare
    ///
    /// range-based version of the lexicographic_compare std algorithm
    ///
    /// \pre Rng1 meets the requirements for an Input range
    /// \pre Rng2 meets the requirements for an Input range
    template<typename Rng1,typename Rng2>
    inline bool lexicographical_compare(Rng1 const & rng1,Rng2 const & rng2)
    {
        return std::lexicographical_compare(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                                            range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2));
    }

    /// \overload
    template<typename Rng1,typename Rng2,typename BinPred>
    inline bool lexicographical_compare(Rng1 const & rng1,Rng2 const & rng2,BinPred pred)
    {
        return std::lexicographical_compare(range_ex_detail::adl_begin(rng1),range_ex_detail::adl_end(rng1),
                                            range_ex_detail::adl_begin(rng2),range_ex_detail::adl_end(rng2),pred);
    }

    /////////////////////////////////////////////////////////////////////////
    // Permutations
    /////////////////////////////////////////////////////////////////////////

    /// \brief template function next_permutation
    ///
    /// range-based version of the next_permutation std algorithm
    ///
    /// \pre Rng meets the requirements for a Bidirectional range
    template<typename Rng>
    inline bool next_permutation(Rng & rng)
    {
        return std::next_permutation(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng>
    inline bool next_permutation(Rng const & rng)
    {
        return std::next_permutation(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng,typename Cmp>
    inline bool next_permutation(Rng & rng,Cmp cmp)
    {
        return std::next_permutation(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),cmp);
    }

    /// \overload
    template<typename Rng,typename Cmp>
    inline bool next_permutation(Rng const & rng,Cmp cmp)
    {
        return std::next_permutation(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),cmp);
    }

    /// \brief template function prev_permutation
    ///
    /// range-based version of the prev_permutation std algorithm
    ///
    /// \pre Rng meets the requirements for a Bidirectional range
    template<typename Rng>
    inline bool prev_permutation(Rng & rng)
    {
        return std::prev_permutation(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng>
    inline bool prev_permutation(Rng const & rng)
    {
        return std::prev_permutation(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng));
    }

    /// \overload
    template<typename Rng,typename Cmp>
    inline bool prev_permutation(Rng & rng,Cmp cmp)
    {
        return std::prev_permutation(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),cmp);
    }

    /// \overload
    template<typename Rng,typename Cmp>
    inline bool prev_permutation(Rng const & rng,Cmp cmp)
    {
        return std::prev_permutation(range_ex_detail::adl_begin(rng),range_ex_detail::adl_end(rng),cmp);
    }

}

#endif
