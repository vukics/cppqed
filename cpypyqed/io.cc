// Copyright Raimar Sandner 2012–2019.
// Copyright András Vukics 2019–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
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
#include "Trajectory.tcc"

#include "blitz2numpy.h"

#include <boost/preprocessor/iteration/local.hpp>

#include <algorithm>

using namespace trajectory;
using namespace boost::python;

namespace cpypyqed { // helps resolving function overloads

template<typename T, int RANK>
object doRead(std::shared_ptr<std::istream> ifs)
{
  using ARRAY=blitz::Array<T,RANK>;
  
  list states, times;
  ARRAY a{ExtTiny<RANK>(0ul)};
  AdaptiveIO<ARRAY> traj{evolved::makeIO(std::move(a))};
  while ( (ifs->peek(), !ifs->eof()) ) {
    trajectory::readViaSStream(traj,ifs);
    states.append(cpputils::arrayToNumpy(a));
    times.append(traj.getTime());
  }
  return make_tuple(states,times);
}

template<typename T, int RANK>
void doWrite(std::shared_ptr<std::ostream> ofs, const numpy::ndarray &a, double time)
{
  trajectory::writeViaSStream(AdaptiveIO<blitz::Array<T,RANK>>{evolved::makeIO(cpputils::numpyToArray<T,RANK>(a), time)},ofs);
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

  trajectory::SerializationMetadata meta{trajectory::readMeta(trajectory::openStateFileReading(f))};

  throw_rank(meta.rank);
  throw_type(meta.typeID);

  std::shared_ptr<std::istream> is = trajectory::openStateFileReading(f);

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
  auto ofs = trajectory::openStateFileWriting(f,std::ios_base::trunc | std::ios_base::binary);

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

  class_<trajectory::SerializationMetadata>("SerializationMetadata")
    .def_readonly("protocolVersion", &trajectory::SerializationMetadata::protocolVersion)
    .def_readonly("rank",            &trajectory::SerializationMetadata::rank)
    .def_readonly("trajectoryID",    &trajectory::SerializationMetadata::trajectoryID);
  
}
