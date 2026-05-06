// Copyright András Vukics 2006–2026. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Trajectory.h"

#include <fstream>


using namespace std;


struct StateFileOpeningException : runtime_error
{
  StateFileOpeningException(const string& filename) : runtime_error("State file opening error: "+filename) {}
};


shared_ptr<istream> cppqedutils::trajectory::openStateFileReading(const string &filename)
{
  shared_ptr<ifstream> ifs = make_shared<ifstream>(filename, ios_base::binary);
  if (!ifs->is_open()) throw StateFileOpeningException(filename);
  return ifs;
}

shared_ptr<ostream> cppqedutils::trajectory::openStateFileWriting(const string &filename, const ios_base::openmode mode)
{
  shared_ptr<ofstream> ofs = make_shared<ofstream>(filename, mode | ios_base::binary);
  if (!ofs->is_open()) throw StateFileOpeningException(filename);
  return ofs;
}
