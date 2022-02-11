// Copyright Raimar Sandner 2012–2019.
// Copyright András Vukics 2019–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef CPYPYQED_BOOST_PYTHON_MODULE_NAME
#ifndef NDEBUG
#define CPYPYQED_BOOST_PYTHON_MODULE_NAME cpypyqed_d
#else // NDEBUG
#define CPYPYQED_BOOST_PYTHON_MODULE_NAME cpypyqed
#endif
#endif // CPYPYQED_BOOST_PYTHON_MODULE_NAME

#define PYTHON_MAX_RANK BOOST_PYTHON_MAX_ARITY

#if PYTHON_MAX_RANK > BLITZ_ARRAY_LARGEST_RANK
#define BLITZ_ARRAY_LARGEST_RANK PYTHON_MAX_RANK
#endif

#include "BlitzArray.h"
#include "Trajectory.h"

#include "blitz2numpy.h"

#include <boost/preprocessor/iteration/local.hpp>

#include <algorithm>

using namespace cppqedutils::trajectory;
using namespace boost::python;


namespace cppqedutils::trajectory {

template <>
struct ReadState<SerializationMetadata>
{
  
  static iarchive& _(SerializationMetadata& sm, iarchive& iar)
  {
    iar & sm;
    return iar;
  }
  
};


template <typename T, int RANK>
using AdaptiveIO = std::tuple<blitz::Array<T,RANK>,double>;


template <typename T, int RANK>
struct ReadState<AdaptiveIO<T,RANK>>
{
  
  static iarchive& _(AdaptiveIO<T,RANK>& arrayTime, iarchive& iar)
  {
    SerializationMetadata sm;
    iar & sm & std::get<0>(arrayTime) & std::get<1>(arrayTime);
    return iar;
  }
  
};


template <typename T, int RANK>
struct WriteState<AdaptiveIO<T,RANK>>
{
  
  static oarchive& _(const AdaptiveIO<T,RANK>& arrayTime, oarchive& oar)
  {
    oar & SerializationMetadata{"CArray",SerializationMetadata::ARRAY_ONLY,RANK} & std::get<0>(arrayTime) & std::get<1>(arrayTime);
    return oar;
  }
  
};

} // cppqedutils::trajectory

namespace cpypyqed { // helps resolving function overloads

template<typename T, int RANK>
object doRead(std::shared_ptr<std::istream> ifs)
{
  list states, times;
  
  AdaptiveIO<T,RANK> arrayTime{};
  
  while ( (ifs->peek(), !ifs->eof()) ) {
    readViaSStream(arrayTime,ifs);
    // specialize cppqedutils::trajectory::readState for something that eats a SerializationMetadata & then reads an array
    // (& ignores the rest of the archive)
    states.append(cppqedutils::arrayToNumpy(std::get<0>(arrayTime)));
    times.append(std::get<1>(arrayTime));
  }
  return make_tuple(states,times);
}

template<typename T, int RANK>
void doWrite(std::shared_ptr<std::ostream> ofs, const numpy::ndarray &a, double time)
{
  AdaptiveIO<T,RANK> arrayTime{cppqedutils::numpyToArray<T,RANK>(a), time};
  writeViaSStream(arrayTime,ofs);
}


void throw_rank(int r)
{
  if (r>PYTHON_MAX_RANK){
    PyErr_SetString(PyExc_NotImplementedError,(std::to_string(r)+">PYTHON_MAX_RANK="+std::to_string(PYTHON_MAX_RANK)).c_str());
    throw_error_already_set();
  }
}


void throw_type(std::string typeID)
{
  if(typeID!="CArray" && typeID!="DArray") {
    PyErr_SetString(PyExc_NotImplementedError, (typeID+": import of this data type not implemented.").c_str());
    throw_error_already_set();
  }
}


object read(str filename) 
{
  auto f = extract<std::string>(filename);

  SerializationMetadata meta;
  readViaSStream(meta,openStateFileReading(f)); // necessary, since this is the only way to determine the rank for the switch statement below
    
  throw_rank(meta.rank);
  throw_type(meta.typeID);

  std::shared_ptr<std::istream> is = openStateFileReading(f);

  list result;
  result.append(meta);

  switch (meta.rank) {
    #define BOOST_PP_LOCAL_MACRO(n) \
    case n: \
      if(meta.typeID=="CArray") result.extend(doRead<dcomp ,n>(is)); \
      if(meta.typeID=="DArray") result.extend(doRead<double,n>(is));  \
      break;
    #define BOOST_PP_LOCAL_LIMITS (1, PYTHON_MAX_RANK)
    #include BOOST_PP_LOCAL_ITERATE()
  }
  return std::move(result);
}

void write(str filename, const numpy::ndarray &array, double time)
{
  auto f = extract<std::string>(filename);
  auto ofs = openStateFileWriting(f,std::ios_base::trunc | std::ios_base::binary);

  const int r=array.get_nd(); throw_rank(r);
  
  switch (r) {
    #define BOOST_PP_LOCAL_MACRO(n)                 \
    case n:                                       \
      if(array.get_dtype()==numpy::dtype::get_builtin<double>())  \
        doWrite<double,n>(ofs,array,time);  \
      if(array.get_dtype()==numpy::dtype::get_builtin<dcomp >()) \
        doWrite<dcomp ,n>(ofs,array,time);   \
      break;
    #define BOOST_PP_LOCAL_LIMITS (1, PYTHON_MAX_RANK)
    #include BOOST_PP_LOCAL_ITERATE()
  }
}

} // cpypyqed


BOOST_PYTHON_MODULE(CPYPYQED_BOOST_PYTHON_MODULE_NAME) {  // Thing in brackets should match output library name
  Py_Initialize();
  numpy::initialize();
  
  def("read", cpypyqed::read,
R"doc(Read in a state vector file.

:param str fname: The input filename.
:returns: A tuple of the form :samp:`(meta, states, times)`.)doc",
  arg("fname")
  );

  def("write", cpypyqed::write,
R"doc(Write a state vector file.

:param str fname: The output filename.
:param ndarray a: The array to write.
:param double t: The time.)doc",
  (arg("fname"),"a", "t")
  );

  class_<SerializationMetadata>("SerializationMetadata")
    .def_readonly("protocolVersion", &SerializationMetadata::protocolVersion)
    .def_readonly("rank",            &SerializationMetadata::rank)
    .def_readonly("trajectoryID",    &SerializationMetadata::trajectoryID);
  
}
