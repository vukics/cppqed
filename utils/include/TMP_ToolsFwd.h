// -*- C++ -*-
#ifndef   UTILS_INCLUDE_TMP_TOOLSFWD_H_INCLUDED
#define   UTILS_INCLUDE_TMP_TOOLSFWD_H_INCLUDED


namespace tmptools {


template<typename, bool IS_CONST>
struct ConditionalAddConst;


template<int N, int Nbeg>
struct Range;

template<int> struct Ordinals;


template<int, int> struct Power;


template<typename Seq, typename ICW> // ICW stands for integral constant wrapper
struct numerical_contains;

template<typename Seq, typename T, T VALUE>
struct numerical_contains_c;


template<int>
struct IsEvenAssert;


template<int, int, bool>
struct pair_c;

} // tmptools


#endif // UTILS_INCLUDE_TMP_TOOLSFWD_H_INCLUDED
