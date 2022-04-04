// Copyright András Vukics 2020–2022.
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

#include <boost/hana.hpp>

#include <cstddef>
#include <iterator>
#include <memory>
#include <iostream>
#include <set>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <valarray>
#include <vector>

#include <algorithm> 
#include <locale>
#include <utility>


// Delimiter definitions with templated variables!

namespace hana=boost::hana;

namespace container_io {

namespace detail {

// Trimming implemented from here: https://stackoverflow.com/a/217605/1171157 in a simplified way
  
static inline auto ltrim(std::string s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
    return !std::isspace(ch);
  }));
  return s;
}

static inline auto rtrim(std::string s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
    return !std::isspace(ch);
  }).base(), s.end());
  return s;
}

static inline auto trim(std::string s) {
  return ltrim(rtrim((s)));
}


static const auto       has_iterator = hana::is_valid([](auto t) -> hana::type<typename decltype(t)::type::      iterator> { });
static const auto has_const_iterator = hana::is_valid([](auto t) -> hana::type<typename decltype(t)::type::const_iterator> { });

static const auto has_begin = hana::is_valid( [](auto t) -> decltype( (void)hana::traits::declval(t).begin() ) { });
static const auto has_end = hana::is_valid([](auto t) -> decltype( (void)hana::traits::declval(t).end() ) { });

static const auto has_cbegin = hana::is_valid([](auto t) -> decltype( (void)hana::traits::declval(t).cbegin() ) { });
static const auto has_cend = hana::is_valid([](auto t) -> decltype( (void)hana::traits::declval(t).cend() ) { });


}  // namespace detail


// Holds the delimiter values as a string type

template <typename String>
struct delimiters_values
{
  using string_type = String;
  String prefix;
  String delimiter;
  String postfix;
};


// Defines the delimiter values for a specific container and character type

template <typename T, typename String>
struct delimiters
{
  using type = delimiters_values<String>;
  static const type values; 
};

// Functor to read/write containers. You can use this directly if you want to specificy a non-default delimiters type.
// The writing logic can be customized by specializing the nested template.

template <bool IS_OUT,
          typename T,
          typename String = ::std::string,
          typename TDelimiters = delimiters<T, String>>
struct container_io_helper
{
  using delimiters_type = TDelimiters;
  
  using traits_type = typename String::traits_type;
  using   char_type = typename String:: value_type;
  
  static const auto constexpr is_out = hana::bool_c<IS_OUT==true>;
  using stream_type = typename decltype(hana::if_(is_out,
                                                  hana::type_c<std::basic_ostream<char_type,traits_type>>,
                                                  hana::type_c<std::basic_istream<char_type,traits_type>>))::type;

  using container_store_type = typename decltype(hana::if_(is_out,hana::traits::add_const(hana::type_c<T>),hana::type_c<T>))::type;

  template<typename ELEM>
  static void atomic(stream_type& s, ELEM& e)
  {
    if constexpr (IS_OUT)
      s<<e;
    else
      s>>e;
  }
  
  template<typename S> // The nested template here is necessary because of the type-erasure mechanism for custom-delimiter specification below
  static void handleDelimiter(stream_type& stream, S delim)
  {
    if constexpr (IS_OUT) {
      if (!delim.empty()) stream << delim;
    }
    else { // There is some liberality with whitespaces
      auto trimmedDelim{detail::trim(delim)};

      for (; std::isspace(stream.peek()); stream.get()) ; // drop whitespaces from the beginning of the stream as well

      for (auto sc : trimmedDelim) {
        char_type c=stream.get();
        if (c != sc) stream.clear(stream_type::badbit);
      }
    }
  }
  
  template <typename U> // The nested template here is necessary because of the type-erasure mechanism for custom-delimiter specification below
  struct io
  {
    static void io_body(typename decltype(hana::if_(is_out,
                                                    hana::traits::add_const(hana::type_c<U>),
                                                    hana::type_c<U>))::type & c, 
                        stream_type & stream)
    {
      auto it = std::begin(c);
      const auto the_end = std::end(c);

      if (it != the_end) {
        for ( ; ; ) {
          atomic(stream,*it);
          
          if (++it == the_end) break;
          
          handleDelimiter(stream,delimiters_type::values.delimiter);
        }
      }
    }
  };

  container_io_helper(container_store_type & container) : container_(container) {}

  inline void operator()(stream_type & stream) const
  {
    handleDelimiter(stream,delimiters_type::values.prefix);

    io<T>::io_body(container_, stream);

    handleDelimiter(stream,delimiters_type::values.postfix);
  }

private:
  container_store_type & container_;
  
};  


// Specialization for pairs

template <bool IS_OUT, typename T, typename String, typename TDelimiters>
template <typename T1, typename T2>
struct container_io_helper<IS_OUT, T, String, TDelimiters>::io<std::pair<T1, T2>>
{
  using container_type_here = typename decltype(hana::if_(is_out,hana::traits::add_const(hana::type_c<std::pair<T1, T2>>),hana::type_c<std::pair<T1, T2>>))::type;
  
  static void io_body(container_type_here & c, stream_type & stream)
  {
    atomic(stream,c.first);
    handleDelimiter(stream,delimiters_type::values.delimiter);
    atomic(stream,c.second);
  }
};

// Specialization for tuples

template <bool IS_OUT, typename T, typename String, typename TDelimiters>
template <typename ...Args>
struct container_io_helper<IS_OUT, T, String, TDelimiters>::io<std::tuple<Args...>>
{
  using container_type_here = typename decltype(hana::if_(is_out,hana::traits::add_const(hana::type_c<std::tuple<Args...>>),hana::type_c<std::tuple<Args...>>))::type;

  template <size_t N>
  static void tuple_io(container_type_here & c, stream_type & stream)
  {
    if constexpr (N!=sizeof...(Args)) {
      if constexpr (N!=0) {handleDelimiter(stream,delimiters_type::values.delimiter);}
      atomic(stream,std::get<N>(c));
      tuple_io<N+1>(c, stream);
    }
  }
  
  static void io_body(container_type_here & c, stream_type & stream)
  {
    tuple_io<0>(c, stream);
  }

};


// Specialization for hana::tuple

template <bool IS_OUT, typename T, typename String, typename TDelimiters>
template <typename ...Args>
struct container_io_helper<IS_OUT, T, String, TDelimiters>::io<hana::tuple<Args...>>
{
  using container_type_here = typename decltype(hana::if_(is_out,hana::traits::add_const(hana::type_c<hana::tuple<Args...>>),hana::type_c<hana::tuple<Args...>>))::type;

  static void io_body(container_type_here & c, stream_type & stream)
  {
    unsigned i=0;
    
    hana::for_each(c, [&](auto& member) {
      if (i) handleDelimiter(stream,delimiters_type::values.delimiter);
      ++i;
      atomic(stream,member);
    });
  }

};


// Prints a write_container_helper to the specified stream.

template<typename T, typename String, typename TDelimiters>
inline auto &
operator<<(std::basic_ostream<typename String::value_type,typename String::traits_type> & stream, const container_io_helper<true, T, String, TDelimiters> & helper)
{
  helper(stream);
  return stream;
}

// reads a read_container_helper from the specified stream.

template<typename T, typename String, typename TDelimiters>
inline auto &
operator>>(std::basic_istream<typename String::value_type,typename String::traits_type> & stream, container_io_helper<false, T, String, TDelimiters> & helper)
{
  helper(stream);
  return stream;
}


// Basic is_container template; specialize to derive from std::true_type for all desired container types

template <typename T>
struct is_container : public std::integral_constant<bool,
                                                    (detail::has_const_iterator(hana::type_c<T>) && detail::has_cbegin(hana::type_c<T>) && detail::has_cend(hana::type_c<T>)) ||
                                                    (detail::has_iterator(hana::type_c<T>) && detail::has_begin(hana::type_c<T>) && detail::has_end(hana::type_c<T>))
                                                    > { };

template <typename T, std::size_t N>
struct is_container<T[N]> : std::true_type { };

template <std::size_t N>
struct is_container<char[N]> : std::false_type { };

template <typename T>
struct is_container<std::valarray<T>> : std::true_type { };

template <typename T1, typename T2>
struct is_container<std::pair<T1, T2>> : std::true_type { };

template <typename ...Args>
struct is_container<std::tuple<Args...>> : std::true_type { };

template <typename ...Args>
struct is_container<hana::tuple<Args...>> : std::true_type { };


// Default delimiters

template <typename T> struct delimiters<T, ::std::string> { static const delimiters_values<::std::string> values; };
template <typename T> const delimiters_values<::std::string> delimiters<T, ::std::string>::values = { "[", ", ", "]" };
template <typename T> struct delimiters<T, ::std::wstring> { static const delimiters_values<::std::wstring> values; };
template <typename T> const delimiters_values<::std::wstring> delimiters<T, ::std::wstring>::values = { L"[", L", ", L"]" };


// Delimiters for (multi)set and unordered_(multi)set

template <typename T, typename TComp, typename TAllocator>
struct delimiters< ::std::set<T, TComp, TAllocator>, ::std::string> { static const delimiters_values<::std::string> values; };

template <typename T, typename TComp, typename TAllocator>
const delimiters_values<::std::string> delimiters< ::std::set<T, TComp, TAllocator>, ::std::string>::values = { "{", ", ", "}" };

template <typename T, typename TComp, typename TAllocator>
struct delimiters< ::std::set<T, TComp, TAllocator>, ::std::wstring> { static const delimiters_values<::std::wstring> values; };

template <typename T, typename TComp, typename TAllocator>
const delimiters_values<::std::wstring> delimiters< ::std::set<T, TComp, TAllocator>, ::std::wstring>::values = { L"{", L", ", L"}" };

template <typename T, typename TComp, typename TAllocator>
struct delimiters< ::std::multiset<T, TComp, TAllocator>, ::std::string> { static const delimiters_values<::std::string> values; };

template <typename T, typename TComp, typename TAllocator>
const delimiters_values<::std::string> delimiters< ::std::multiset<T, TComp, TAllocator>, ::std::string>::values = { "{", ", ", "}" };

template <typename T, typename TComp, typename TAllocator>
struct delimiters< ::std::multiset<T, TComp, TAllocator>, ::std::wstring> { static const delimiters_values<::std::wstring> values; };

template <typename T, typename TComp, typename TAllocator>
const delimiters_values<::std::wstring> delimiters< ::std::multiset<T, TComp, TAllocator>, ::std::wstring>::values = { L"{", L", ", L"}" };

template <typename T, typename THash, typename TEqual, typename TAllocator>
struct delimiters< ::std::unordered_set<T, THash, TEqual, TAllocator>, ::std::string> { static const delimiters_values<::std::string> values; };

template <typename T, typename THash, typename TEqual, typename TAllocator>
const delimiters_values<::std::string> delimiters< ::std::unordered_set<T, THash, TEqual, TAllocator>, ::std::string>::values = { "{", ", ", "}" };

template <typename T, typename THash, typename TEqual, typename TAllocator>
struct delimiters< ::std::unordered_set<T, THash, TEqual, TAllocator>, ::std::wstring> { static const delimiters_values<::std::wstring> values; };

template <typename T, typename THash, typename TEqual, typename TAllocator>
const delimiters_values<::std::wstring> delimiters< ::std::unordered_set<T, THash, TEqual, TAllocator>, ::std::wstring>::values = { L"{", L", ", L"}" };

template <typename T, typename THash, typename TEqual, typename TAllocator>
struct delimiters< ::std::unordered_multiset<T, THash, TEqual, TAllocator>, ::std::string> { static const delimiters_values<::std::string> values; };

template <typename T, typename THash, typename TEqual, typename TAllocator>
const delimiters_values<::std::string> delimiters< ::std::unordered_multiset<T, THash, TEqual, TAllocator>, ::std::string>::values = { "{", ", ", "}" };

template <typename T, typename THash, typename TEqual, typename TAllocator>
struct delimiters< ::std::unordered_multiset<T, THash, TEqual, TAllocator>, ::std::wstring> { static const delimiters_values<::std::wstring> values; };

template <typename T, typename THash, typename TEqual, typename TAllocator>
const delimiters_values<::std::wstring> delimiters< ::std::unordered_multiset<T, THash, TEqual, TAllocator>, ::std::wstring>::values = { L"{", L", ", L"}" };


// Delimiters for pair and tuple

template <typename T1, typename T2> struct delimiters<std::pair<T1, T2>, ::std::string> { static const delimiters_values<::std::string> values; };
template <typename T1, typename T2> const delimiters_values<::std::string> delimiters<std::pair<T1, T2>, ::std::string>::values = { "(", ", ", ")" };
template <typename T1, typename T2> struct delimiters< ::std::pair<T1, T2>, ::std::wstring> { static const delimiters_values<::std::wstring> values; };
template <typename T1, typename T2> const delimiters_values<::std::wstring> delimiters< ::std::pair<T1, T2>, ::std::wstring>::values = { L"(", L", ", L")" };

template <typename ...Args> struct delimiters<std::tuple<Args...>, ::std::string> { static const delimiters_values<::std::string> values; };
template <typename ...Args> const delimiters_values<::std::string> delimiters<std::tuple<Args...>, ::std::string>::values = { "(", ", ", ")" };
template <typename ...Args> struct delimiters< ::std::tuple<Args...>, ::std::wstring> { static const delimiters_values<::std::wstring> values; };
template <typename ...Args> const delimiters_values<::std::wstring> delimiters< ::std::tuple<Args...>, ::std::wstring>::values = { L"(", L", ", L")" };


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
    return s << container_io_helper<true, T, ::std::string, Delims>(t);
  }

  std::wostream & stream(std::wostream & s)
  {
    return s << container_io_helper<true, T, ::std::wstring, Delims>(t);
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
inline auto & operator<<(std::basic_ostream<TChar, TCharTraits> & s, const custom_delims<Delims> & p)
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
  inline const_iterator cbegin() const { return _array; }
  inline const_iterator cend() const { return _array + _n; }
  inline const_iterator begin() const { return cbegin(); }
  inline const_iterator end() const { return cend(); }

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

  const_iterator cbegin() const
  {
    return m_map.cbegin(n);
  }

  const_iterator cend() const
  {
    return m_map.cend(n);
  }

  inline const_iterator begin() const { return cbegin(); }
  inline const_iterator end() const { return cend(); }
  
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
  return stream << ::container_io::container_io_helper<true,T, basic_string<TChar,TCharTraits>>(container);
}

template<typename T, typename TChar, typename TCharTraits>
inline typename enable_if< ::container_io::is_container<T>::value, basic_istream<TChar, TCharTraits> &>::type
operator>>(basic_istream<TChar, TCharTraits> & stream, T & container)
{
  auto c{::container_io::container_io_helper<false, T, basic_string<TChar,TCharTraits>>(container)};
  return stream >> c;
}

} // std


namespace boost { namespace hana
{

// Prints a container to the stream using default delimiters

template<typename TChar, typename TCharTraits, typename ...Args>
inline auto &
operator<<(::std::basic_ostream<TChar, TCharTraits> & stream, const hana::tuple<Args...> & container)
{
  return stream << ::container_io::container_io_helper<true, hana::tuple<Args...>, ::std::basic_string<TChar,TCharTraits>>(container);
}

template<typename TChar, typename TCharTraits, typename ...Args>
inline auto &
operator>>(::std::basic_istream<TChar, TCharTraits> & stream, hana::tuple<Args...> & container)
{
  auto c{::container_io::container_io_helper<false, hana::tuple<Args...>, ::std::basic_string<TChar,TCharTraits>>(container)};
  return stream >> c;
}

} } // hana


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
