// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ArrayTraits.h"
#include "BlitzArray.h"
#include "Trajectory.tcc"

#ifndef DO_NOT_USE_BOOST_COMPRESSION
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/device/file.hpp>
#endif // DO_NOT_USE_BOOST_COMPRESSION


using namespace std;

#ifndef    DO_NOT_USE_BOOST_COMPRESSION
namespace {

bool isbz2(const std::string filename)
{
  using namespace std;
  ifstream file(filename, ios_base::in | ios_base::binary);
  if (file.peek() == ifstream::traits_type::eof())
    return true;
  string header; header.resize(3);
  file.read(header.data(),3);
  return header=="BZh";
}

}

#endif // DO_NOT_USE_BOOST_COMPRESSION


struct StateFileOpeningException : std::runtime_error
{
  StateFileOpeningException(const std::string& filename) : std::runtime_error("State file opening error: "+filename) {}
};


std::shared_ptr<istream> trajectory::openStateFileReading(const std::string &filename)
{
#ifdef DO_NOT_USE_BOOST_COMPRESSION
  std::shared_ptr<ifstream> ifs = std::make_shared<ifstream>(filename, ios_base::in | ios_base::binary);
  if (!ifs->is_open()) throw StateFileOpeningException(filename);
  return ifs;
#else
  using namespace boost::iostreams;
  std::shared_ptr<filtering_istream> in = std::make_shared<filtering_istream>();
  file_source file(filename, std::ios_base::in | std::ios_base::binary);
  if (!file.is_open()) throw StateFileOpeningException(filename);
  if (isbz2(filename)) in->push(bzip2_decompressor());
  in->push(file);
  return in;
#endif // DO_NOT_USE_BOOST_COMPRESSION
}

std::shared_ptr<ostream> trajectory::openStateFileWriting(const std::string &filename, const ios_base::openmode mode)
{
#ifdef DO_NOT_USE_BOOST_COMPRESSION
  std::shared_ptr<ofstream> ofs = std::make_shared<ofstream>(filename, mode | std::ios_base::binary );
  if (!ofs->is_open()) throw StateFileOpeningException(filename);
  return ofs;
#else
  using namespace boost::iostreams;
  std::shared_ptr<filtering_ostream> out = std::make_shared<filtering_ostream>();
  file_sink file(filename, mode | std::ios_base::binary );
  if (!file.is_open()) throw StateFileOpeningException(filename);
  if (isbz2(filename)) out->push(bzip2_compressor());
  out->push(file);
  return out;
#endif // DO_NOT_USE_BOOST_COMPRESSION
}


void trajectory::details::writeNextArchive(ostream* ofs, const ostringstream &oss)
{
  const string& buffer=oss.str();
  *ofs<<buffer.size(); ofs->write(&buffer[0],buffer.size());
}

void trajectory::details::readNextArchive(istream* ifs, istringstream &iss)
{
  string buffer;
  streamsize n; *ifs>>n; buffer.resize(n);
  ifs->read(&buffer[0],n);
  iss.str(buffer);
}


auto trajectory::readMeta(istream* ifs) -> SerializationMetadata
{
  istringstream iss(ios_base::binary);
  details::readNextArchive(ifs,iss);
  cpputils::iarchive archive(iss);
  SerializationMetadata meta;
  archive >> meta;
  return meta;
}

const std::string trajectory::SerializationMetadata::UNSPECIFIED = "Unspecified";
const std::string trajectory::SerializationMetadata::ARRAY_ONLY  = "ArrayOnly";

