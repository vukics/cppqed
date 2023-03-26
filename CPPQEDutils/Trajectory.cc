// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Trajectory.h"

#ifndef DO_NOT_USE_BOOST_COMPRESSION
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/device/file.hpp>
#endif // DO_NOT_USE_BOOST_COMPRESSION


using namespace std;

#ifndef    DO_NOT_USE_BOOST_COMPRESSION
namespace {

bool isbz2(const string filename)
{
  using namespace std;
  ifstream file(filename, ios_base::binary);
  if (file.peek() == ifstream::traits_type::eof())
    return true;
  string header; header.resize(3);
  file.read(header.data(),3);
  return header=="BZh";
}

}

#endif // DO_NOT_USE_BOOST_COMPRESSION


struct StateFileOpeningException : runtime_error
{
  StateFileOpeningException(const string& filename) : runtime_error("State file opening error: "+filename) {}
};


shared_ptr<istream> cppqedutils::trajectory::openStateFileReading(const string &filename)
{
#ifdef DO_NOT_USE_BOOST_COMPRESSION
  shared_ptr<ifstream> ifs = make_shared<ifstream>(filename, ios_base::binary);
  if (!ifs->is_open()) throw StateFileOpeningException(filename);
  return ifs;
#else
  using namespace boost::iostreams;
  shared_ptr<filtering_istream> in = make_shared<filtering_istream>();
  file_source file(filename, ios_base::in | ios_base::binary);
  if (!file.is_open()) throw StateFileOpeningException(filename);
  if (isbz2(filename)) in->push(bzip2_decompressor());
  in->push(file);
  return in;
#endif // DO_NOT_USE_BOOST_COMPRESSION
}

shared_ptr<ostream> cppqedutils::trajectory::openStateFileWriting(const string &filename, const ios_base::openmode mode)
{
#ifdef DO_NOT_USE_BOOST_COMPRESSION
  shared_ptr<ofstream> ofs = make_shared<ofstream>(filename, mode | ios_base::binary );
  if (!ofs->is_open()) throw StateFileOpeningException(filename);
  return ofs;
#else
  using namespace boost::iostreams;
  shared_ptr<filtering_ostream> out = make_shared<filtering_ostream>();
  file_sink file(filename, mode | ios_base::binary );
  if (!file.is_open()) throw StateFileOpeningException(filename);
  if (isbz2(filename)) out->push(bzip2_compressor());
  out->push(file);
  return out;
#endif // DO_NOT_USE_BOOST_COMPRESSION
}


std::ostream& cppqedutils::trajectory::stream(double tdp, std::ostream& os) {return os<<tdp;}
std::ostream& cppqedutils::trajectory::stream(dcomp  tdp, std::ostream& os) {return os<<tdp;}



