// Copyright András Vukics 2020
// Extending the original work of Louis Delacroix with operator>> + boost::hana::tuple
// Code retrieved from https://github.com/louisdx/cxx-prettyprint on 13/03/2020

//          Copyright Louis Delacroix 2010 - 2014.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)
//
// A pretty printing library for C++
//
// Usage:
// Include this header, and operator<< will "just work".


#ifndef CONTAINER_IO_H_INCLUDED
#define CONTAINER_IO_H_INCLUDED

#include <boost/hana/tuple.hpp>
#include <boost/hana/for_each.hpp>

#include <cstddef>
#include <iterator>
#include <memory>
#include <ostream>
#include <set>
#include <tuple>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <valarray>
#include <vector>

namespace container_io {

namespace detail {
// SFINAE type trait to detect whether T::const_iterator exists.

struct sfinae_base
{
  using yes = char;
  using no  = yes[2];
};

template <typename T>
struct has_const_iterator : private sfinae_base
{
private:
  template <typename C> static yes & test(typename C::const_iterator*);
  template <typename C> static no  & test(...);
public:
  static const bool value = sizeof(test<T>(nullptr)) == sizeof(yes);
  using type =  T;
};

template <typename T>
struct has_begin_end : private sfinae_base
{
private:
  template <typename C>
  static yes & f(typename std::enable_if<std::is_same<decltype(static_cast<typename C::const_iterator(C::*)() const>(&C::begin)),typename C::const_iterator(C::*)() const>::value>::type *);

  template <typename C> static no & f(...);

  template <typename C>
  static yes & g(typename std::enable_if<std::is_same<decltype(static_cast<typename C::const_iterator(C::*)() const>(&C::end)),typename C::const_iterator(C::*)() const>::value, void>::type*);

  template <typename C> static no & g(...);

public:
  static bool const beg_value = sizeof(f<T>(nullptr)) == sizeof(yes);
  static bool const end_value = sizeof(g<T>(nullptr)) == sizeof(yes);
};

}  // namespace detail


// Holds the delimiter values for a specific character type

template <typename TChar>
struct delimiters_values
{
  using char_type = TChar;
  const char_type * prefix;
  const char_type * delimiter;
  const char_type * postfix;
};


// Defines the delimiter values for a specific container and character type

template <typename T, typename TChar>
struct delimiters
{
  using type = delimiters_values<TChar>;
  static const type values; 
};


// Functor to write containers. You can use this directly if you want
// to specificy a non-default delimiters type. The writing logic can
// be customized by specializing the nested template.

template <typename T,
          typename TChar = char,
          typename TCharTraits = ::std::char_traits<TChar>,
          typename TDelimiters = delimiters<T, TChar>>
struct write_container_helper
{
  using delimiters_type = TDelimiters;
  using ostream_type = std::basic_ostream<TChar, TCharTraits>;

  template <typename U>
  struct writer
  {
    static void write_body(const U & c, ostream_type & stream)
    {
      auto it = std::begin(c);
      const auto the_end = std::end(c);

      if (it != the_end) {
        for ( ; ; ) {
          stream << *it;
          
          if (++it == the_end) break;
          
          if (delimiters_type::values.delimiter != NULL) stream << delimiters_type::values.delimiter;
        }
      }
    }
  };

  write_container_helper(const T & container) : container_(container) {}

  inline void operator()(ostream_type & stream) const
  {
    if (delimiters_type::values.prefix != NULL) stream << delimiters_type::values.prefix;

    writer<T>::write_body(container_, stream);

    if (delimiters_type::values.postfix != NULL) stream << delimiters_type::values.postfix;
  }

private:
  const T & container_;
};


// Functor to read containers.

template <typename T,
          typename TChar = char,
          typename TCharTraits = ::std::char_traits<TChar>,
          typename TDelimiters = delimiters<T, TChar>>
struct read_container_helper
{
  using delimiters_type = TDelimiters;
  using istream_type = std::basic_istream<TChar, TCharTraits>;

  template <typename U>
  struct reader
  {
    static void read_body(U & c, istream_type & stream)
    {
      auto it = std::begin(c);
      const auto the_end = std::end(c);

      if (it != the_end) {
        for ( ; ; ) {
          stream >> *it;
          
          if (++it == the_end) break;
          
          if (delimiters_type::values.delimiter != NULL) { // try to eat delimiter
            TChar c;
            stream >> c;
            if (c!=delimiters_type::values.delimiter) stream.clear(istream_type::badbit);
          }
        }
      }
    }
  };

  read_container_helper(T & container) : container_(container) {}

private:
  void callOperatorHelper(istream_type & stream, const TChar* delim) const
  {
    if (delim != NULL) {
      TChar c;
      stream >> c;
      if (c!=delim) stream.clear(istream_type::badbit);
    }
  }
  
public:
  
  void operator()(istream_type & stream)
  {
    callOperatorHelper(stream,delimiters_type::values.prefix);

    reader<T>::read_body(container_, stream);

    callOperatorHelper(stream,delimiters_type::values.postfix);
  }

private:
  // What kind of reference to store here?
  T & container_;
};


// Specialization for pairs

template <typename T, typename TChar, typename TCharTraits, typename TDelimiters>
template <typename T1, typename T2>
struct write_container_helper<T, TChar, TCharTraits, TDelimiters>::writer<std::pair<T1, T2>>
{
  using ostream_type = typename write_container_helper<T, TChar, TCharTraits, TDelimiters>::ostream_type;

  static void write_body(const std::pair<T1, T2> & c, ostream_type & stream)
  {
    stream << c.first;
    if (write_container_helper<T, TChar, TCharTraits, TDelimiters>::delimiters_type::values.delimiter != NULL)
        stream << write_container_helper<T, TChar, TCharTraits, TDelimiters>::delimiters_type::values.delimiter;
    stream << c.second;
  }
};

// Specialization for tuples

template <typename T, typename TChar, typename TCharTraits, typename TDelimiters>
template <typename ...Args>
struct write_container_helper<T, TChar, TCharTraits, TDelimiters>::writer<std::tuple<Args...>>
{
  using ostream_type = typename write_container_helper<T, TChar, TCharTraits, TDelimiters>::ostream_type;
  using element_type = std::tuple<Args...>;

  template <std::size_t I> struct Int { };

  static void write_body(const element_type & c, ostream_type & stream)
  {
    tuple_write(c, stream, Int<0>());
  }

  static void tuple_write(const element_type &, ostream_type &, Int<sizeof...(Args)>) {}

  static void tuple_write(const element_type & c, ostream_type & stream,
                          typename std::conditional<sizeof...(Args) != 0, Int<0>, std::nullptr_t>::type)
  {
    stream << std::get<0>(c);
    tuple_write(c, stream, Int<1>());
  }

  template <std::size_t N>
  static void tuple_write(const element_type & c, ostream_type & stream, Int<N>)
  {
    if (write_container_helper<T, TChar, TCharTraits, TDelimiters>::delimiters_type::values.delimiter != NULL)
      stream << write_container_helper<T, TChar, TCharTraits, TDelimiters>::delimiters_type::values.delimiter;

    stream << std::get<N>(c);

    tuple_write(c, stream, Int<N + 1>());
  }
};


// Specialization for hana::tuple

template <typename T, typename TChar, typename TCharTraits, typename TDelimiters>
template <typename ...Args>
struct write_container_helper<T, TChar, TCharTraits, TDelimiters>::writer<boost::hana::tuple<Args...>>
{
  using ostream_type = typename write_container_helper<T, TChar, TCharTraits, TDelimiters>::ostream_type;
  using element_type = boost::hana::tuple<Args...>;

  static void write_body(const element_type & c, ostream_type & stream)
  {
    unsigned i=0;
    
    boost::hana::for_each(c, [&](const auto& member) {
      if (i && write_container_helper<T, TChar, TCharTraits, TDelimiters>::delimiters_type::values.delimiter != NULL) stream << write_container_helper<T, TChar, TCharTraits, TDelimiters>::delimiters_type::values.delimiter;
      ++i;
      stream << member ;
    });
  }

};


// Prints a write_container_helper to the specified stream.

template<typename T, typename TChar, typename TCharTraits, typename TDelimiters>
inline std::basic_ostream<TChar, TCharTraits> & operator<<(std::basic_ostream<TChar, TCharTraits> & stream, const write_container_helper<T, TChar, TCharTraits, TDelimiters> & helper)
{
  helper(stream);
  return stream;
}

// reads a read_container_helper from the specified stream.

template<typename T, typename TChar, typename TCharTraits, typename TDelimiters>
inline std::basic_istream<TChar, TCharTraits> & operator>>(std::basic_istream<TChar, TCharTraits> & stream, read_container_helper<T, TChar, TCharTraits, TDelimiters> & helper)
{
  helper(stream);
  return stream;
}


// Basic is_container template; specialize to derive from std::true_type for all desired container types

template <typename T>
struct is_container : public std::integral_constant<bool,
                                                    detail::has_const_iterator<T>::value &&
                                                    detail::has_begin_end<T>::beg_value  &&
                                                    detail::has_begin_end<T>::end_value> { };

template <typename T, std::size_t N>
struct is_container<T[N]> : std::true_type { };

template <std::size_t N>
struct is_container<char[N]> : std::false_type { };

template <typename T>
struct is_container<std::vector<T>> : std::true_type { };
// TODO: this of course should not be necessary, but the sfinae-based solution above doesn’t seem to work in the non-const container case (which is the case when reading from istream)

template <typename T>
struct is_container<std::valarray<T>> : std::true_type { };

template <typename T1, typename T2>
struct is_container<std::pair<T1, T2>> : std::true_type { };

template <typename ...Args>
struct is_container<std::tuple<Args...>> : std::true_type { };

template <typename ...Args>
struct is_container<boost::hana::tuple<Args...>> : std::true_type { };


// Default delimiters

template <typename T> struct delimiters<T, char> { static const delimiters_values<char> values; };
template <typename T> const delimiters_values<char> delimiters<T, char>::values = { "[", ", ", "]" };
template <typename T> struct delimiters<T, wchar_t> { static const delimiters_values<wchar_t> values; };
template <typename T> const delimiters_values<wchar_t> delimiters<T, wchar_t>::values = { L"[", L", ", L"]" };


// Delimiters for (multi)set and unordered_(multi)set

template <typename T, typename TComp, typename TAllocator>
struct delimiters< ::std::set<T, TComp, TAllocator>, char> { static const delimiters_values<char> values; };

template <typename T, typename TComp, typename TAllocator>
const delimiters_values<char> delimiters< ::std::set<T, TComp, TAllocator>, char>::values = { "{", ", ", "}" };

template <typename T, typename TComp, typename TAllocator>
struct delimiters< ::std::set<T, TComp, TAllocator>, wchar_t> { static const delimiters_values<wchar_t> values; };

template <typename T, typename TComp, typename TAllocator>
const delimiters_values<wchar_t> delimiters< ::std::set<T, TComp, TAllocator>, wchar_t>::values = { L"{", L", ", L"}" };

template <typename T, typename TComp, typename TAllocator>
struct delimiters< ::std::multiset<T, TComp, TAllocator>, char> { static const delimiters_values<char> values; };

template <typename T, typename TComp, typename TAllocator>
const delimiters_values<char> delimiters< ::std::multiset<T, TComp, TAllocator>, char>::values = { "{", ", ", "}" };

template <typename T, typename TComp, typename TAllocator>
struct delimiters< ::std::multiset<T, TComp, TAllocator>, wchar_t> { static const delimiters_values<wchar_t> values; };

template <typename T, typename TComp, typename TAllocator>
const delimiters_values<wchar_t> delimiters< ::std::multiset<T, TComp, TAllocator>, wchar_t>::values = { L"{", L", ", L"}" };

template <typename T, typename THash, typename TEqual, typename TAllocator>
struct delimiters< ::std::unordered_set<T, THash, TEqual, TAllocator>, char> { static const delimiters_values<char> values; };

template <typename T, typename THash, typename TEqual, typename TAllocator>
const delimiters_values<char> delimiters< ::std::unordered_set<T, THash, TEqual, TAllocator>, char>::values = { "{", ", ", "}" };

template <typename T, typename THash, typename TEqual, typename TAllocator>
struct delimiters< ::std::unordered_set<T, THash, TEqual, TAllocator>, wchar_t> { static const delimiters_values<wchar_t> values; };

template <typename T, typename THash, typename TEqual, typename TAllocator>
const delimiters_values<wchar_t> delimiters< ::std::unordered_set<T, THash, TEqual, TAllocator>, wchar_t>::values = { L"{", L", ", L"}" };

template <typename T, typename THash, typename TEqual, typename TAllocator>
struct delimiters< ::std::unordered_multiset<T, THash, TEqual, TAllocator>, char> { static const delimiters_values<char> values; };

template <typename T, typename THash, typename TEqual, typename TAllocator>
const delimiters_values<char> delimiters< ::std::unordered_multiset<T, THash, TEqual, TAllocator>, char>::values = { "{", ", ", "}" };

template <typename T, typename THash, typename TEqual, typename TAllocator>
struct delimiters< ::std::unordered_multiset<T, THash, TEqual, TAllocator>, wchar_t> { static const delimiters_values<wchar_t> values; };

template <typename T, typename THash, typename TEqual, typename TAllocator>
const delimiters_values<wchar_t> delimiters< ::std::unordered_multiset<T, THash, TEqual, TAllocator>, wchar_t>::values = { L"{", L", ", L"}" };


// Delimiters for pair and tuple

template <typename T1, typename T2> struct delimiters<std::pair<T1, T2>, char> { static const delimiters_values<char> values; };
template <typename T1, typename T2> const delimiters_values<char> delimiters<std::pair<T1, T2>, char>::values = { "(", ", ", ")" };
template <typename T1, typename T2> struct delimiters< ::std::pair<T1, T2>, wchar_t> { static const delimiters_values<wchar_t> values; };
template <typename T1, typename T2> const delimiters_values<wchar_t> delimiters< ::std::pair<T1, T2>, wchar_t>::values = { L"(", L", ", L")" };

template <typename ...Args> struct delimiters<std::tuple<Args...>, char> { static const delimiters_values<char> values; };
template <typename ...Args> const delimiters_values<char> delimiters<std::tuple<Args...>, char>::values = { "(", ", ", ")" };
template <typename ...Args> struct delimiters< ::std::tuple<Args...>, wchar_t> { static const delimiters_values<wchar_t> values; };
template <typename ...Args> const delimiters_values<wchar_t> delimiters< ::std::tuple<Args...>, wchar_t>::values = { L"(", L", ", L")" };


// Type-erasing helper class for easy use of custom delimiters.
// Requires TCharTraits = std::char_traits<TChar> and TChar = char or wchar_t, and MyDelims needs to be defined for TChar.
// Usage: "cout << container_io::custom_delims<MyDelims>(x)".

struct custom_delims_base
{
  virtual ~custom_delims_base() { }
  virtual std::ostream & stream(::std::ostream &) = 0;
  virtual std::wostream & stream(::std::wostream &) = 0;
};

template <typename T, typename Delims>
struct custom_delims_wrapper : custom_delims_base
{
  custom_delims_wrapper(const T & t_) : t(t_) { }

  std::ostream & stream(std::ostream & s)
  {
    return s << write_container_helper<T, char, std::char_traits<char>, Delims>(t);
  }

  std::wostream & stream(std::wostream & s)
  {
    return s << write_container_helper<T, wchar_t, std::char_traits<wchar_t>, Delims>(t);
  }

private:
  const T & t;
};

template <typename Delims>
struct custom_delims
{
  template <typename Container>
  custom_delims(const Container & c) : base(new custom_delims_wrapper<Container, Delims>(c)) { }

  std::unique_ptr<custom_delims_base> base;
};

template <typename TChar, typename TCharTraits, typename Delims>
inline std::basic_ostream<TChar, TCharTraits> & operator<<(std::basic_ostream<TChar, TCharTraits> & s, const custom_delims<Delims> & p)
{
  return p.base->stream(s);
}


// A wrapper for a C-style array given as pointer-plus-size.
// Usage: std::cout << write_array(arr, n) << std::endl;

template<typename T>
struct array_wrapper_n
{
  typedef const T * const_iterator;
  typedef T value_type;

  array_wrapper_n(const T * const a, size_t n) : _array(a), _n(n) { }
  inline const_iterator begin() const { return _array; }
  inline const_iterator end() const { return _array + _n; }

private:
  const T * const _array;
  size_t _n;
};


// A wrapper for hash-table based containers that offer local iterators to each bucket.
// Usage: std::cout << bucket_write(m, 4) << std::endl;  (Prints bucket 5 of container m.)

template <typename T>
struct bucket_write_wrapper
{
  typedef typename T::const_local_iterator const_iterator;
  typedef typename T::size_type size_type;

  const_iterator begin() const
  {
    return m_map.cbegin(n);
  }

  const_iterator end() const
  {
    return m_map.cend(n);
  }

  bucket_write_wrapper(const T & m, size_type bucket) : m_map(m), n(bucket) { }

private:
  const T & m_map;
  const size_type n;
};

} // container_io


// Global accessor functions for the convenience wrappers

template<typename T>
inline container_io::array_wrapper_n<T> write_array(const T * const a, size_t n)
{
  return container_io::array_wrapper_n<T>(a, n);
}

template <typename T> container_io::bucket_write_wrapper<T>
bucket_write(const T & m, typename T::size_type n)
{
  return container_io::bucket_write_wrapper<T>(m, n);
}


// Main magic entry point: An overload snuck into namespace std.
// Can we do better?

namespace std
{
// Prints a container to the stream using default delimiters

template<typename T, typename TChar, typename TCharTraits>
inline typename enable_if< ::container_io::is_container<T>::value, basic_ostream<TChar, TCharTraits> &>::type
operator<<(basic_ostream<TChar, TCharTraits> & stream, const T & container)
{
  return stream << ::container_io::write_container_helper<T, TChar, TCharTraits>(container);
}

template<typename T, typename TChar, typename TCharTraits>
inline typename enable_if< ::container_io::is_container<T>::value, basic_istream<TChar, TCharTraits> &>::type
operator>>(basic_istream<TChar, TCharTraits> & stream, T & container)
{
  auto c{::container_io::read_container_helper<T, TChar, TCharTraits>(container)};
  return stream >> c;
}

} // std


namespace boost { namespace hana
{

// Prints a container to the stream using default delimiters

template<typename TChar, typename TCharTraits, typename ...Args>
inline auto &
operator<<(::std::basic_ostream<TChar, TCharTraits> & stream, const boost::hana::tuple<Args...> & container)
{
  return stream << ::container_io::write_container_helper<boost::hana::tuple<Args...>, TChar, TCharTraits>(container);
}

} } // boost::hana


/*
Copyright 2011 Christopher Allen Ogden. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list
      of conditions and the following disclaimer in the documentation and/or other materials
      provided with the distribution.

THIS SOFTWARE IS PROVIDED BY CHRISTOPHER ALLEN OGDEN ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL CHRISTOPHER ALLEN OGDEN OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those of the
authors and should not be interpreted as representing official policies, either expressed
or implied, of Christopher Allen Ogden.
*/


namespace boost {
namespace serialization {

template<uint N>
struct Serialize
{
    template<class Archive, typename... Args>
    static void serialize(Archive & ar, std::tuple<Args...> & t, const unsigned int version)
    {
        ar & std::get<N-1>(t);
        Serialize<N-1>::serialize(ar, t, version);
    }
};

template<>
struct Serialize<0>
{
    template<class Archive, typename... Args>
    static void serialize(Archive & ar, std::tuple<Args...> & t, const unsigned int version)
    {
        (void) ar;
        (void) t;
        (void) version;
    }
};

template<class Archive, typename... Args>
void serialize(Archive & ar, std::tuple<Args...> & t, const unsigned int version)
{
    Serialize<sizeof...(Args)>::serialize(ar, t, version);
}

}
}

#endif  // CONTAINER_IO_H_INCLUDED
