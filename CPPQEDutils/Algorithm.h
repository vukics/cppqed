// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Generic algorithms not found in either STL or Boost}
#ifndef   CPPQEDCORE_UTILS_ALGORITHM_H_INCLUDED
#define   CPPQEDCORE_UTILS_ALGORITHM_H_INCLUDED


namespace cppqedutils {


/// Fills a container by output iterator with concatenated values taken subsequently from the input sequences. \note Concatenation can be expressed as accumulation
template<typename SeqOfSeqs, typename Out_Iterator>
const Out_Iterator
concatenateViaIterator(const SeqOfSeqs& sOs, ///<[in] the sequence containing the input sequences
                       Out_Iterator out)
{
  return std::accumulate(sOs.begin(),sOs.end(),out, [](auto iter, const auto& t){return std::copy(t.begin(),t.end(),iter);});
}


/// Fills a container of the necessary size with concatenated values taken subsequently from the input sequences
template<typename SeqOfSeqs, typename Out>
const Out&
concatenate(const SeqOfSeqs& sOs, ///<[in] the sequence containing the input sequences
            Out& out)
{
  concatenateViaIterator(sOs,out.begin());
  return out;
}


/// Fills an *empty* (default-constructed) container with concatenated values taken subsequently from the input sequences
template<typename Out, typename SeqOfSeqs>
const Out
concatenateGrow(const SeqOfSeqs& sOs)
{
  Out empty;
  concatenateViaIterator(sOs,std::back_inserter(empty));
  return empty;
}


} // cppqedutils


#endif // CPPQEDCORE_UTILS_ALGORITHM_H_INCLUDED
