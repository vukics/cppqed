// -*- C++ -*-
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "Evolved.h"
#include "PythonExtension.h"

#include "Types.h"

#include "BlitzArrayTraits.h"
#include "BlitzTiny.h"
#include "Trajectory.tcc"

#if PYTHON_MAX_RANK > BLITZ_ARRAY_LARGEST_RANK
#define BLITZ_ARRAY_LARGEST_RANK PYTHON_MAX_RANK
#endif

#include <algorithm>
#include <blitz/tinyvec2.h>
#include <boost/function.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/preprocessor/iteration/local.hpp>
#include <boost/python/exception_translator.hpp>
#include <numpy/ndarrayobject.h>
#include <fstream>
#include <iostream>
#include <string>

using namespace trajectory;
using namespace boost::python;

namespace pythonext {

namespace {

template<typename A, int RANK>
object doRead(std::ifstream &ifs, NPY_TYPES npy_dtype)
{
  list states;
  list times;
  list result;
  ExtTiny<RANK> dims;
  dims=0;
  A a(dims);
  AdaptiveIO<A> traj(evolved::makeIO(a));
  while ( (ifs.peek(), !ifs.eof()) ) {
    trajectory::readViaSStream(traj,ifs);
    npy_intp npy_dims[RANK];
    dims = a.extent();
    std::copy(dims.begin(),dims.end(), npy_dims);
    PyObject * pyObj = PyArray_SimpleNewFromData(RANK, npy_dims, npy_dtype, a.dataFirst());
    handle<> h( pyObj );
    numeric::array arr( h );
    states.append(arr.copy());
    times.append(traj.getEvolvedIO()->getTime());
  }
  return make_tuple(states,times);
}

template<typename A, typename dtype, int RANK>
void doWrite(std::ofstream *ofs, PyArrayObject *a, double time)
{
  npy_intp *dims=PyArray_DIMS(a);
  blitz::TinyVector<int,RANK> shape;
  for (int i=0; i<RANK; i++) shape[i]=dims[i];
  A  blitz_a = A(static_cast<dtype *>(PyArray_DATA(a)), shape, blitz::duplicateData);
  AdaptiveIO<A> traj(evolved::makeIO(blitz_a, time));
  trajectory::writeViaSStream(traj,ofs);
}

template<typename T>
void throw_file(const T &s, const std::string &f)
{
  if (!s.is_open()){
    PyErr_SetString(PyExc_IOError, (std::string("Could not open ")+f).c_str());
    throw_error_already_set();
  }
}

void throw_rank(int r)
{
  if (r>PYTHON_MAX_RANK){
    PyErr_SetString(PyExc_NotImplementedError, (boost::lexical_cast<std::string>(r)+std::string(">PYTHON_MAX_RANK.")).c_str());
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

}

object read(str filename) 
{
  std::string f = extract<std::string>(filename);

  std::ifstream ifs(f.c_str(),std::ios_base::binary);
  throw_file(ifs,f);

  trajectory::SerializationMetadata meta = trajectory::readMeta(ifs);
  
  throw_rank(meta.rank);
  throw_type(meta.typeID);

  list result;
  result.append(meta);

  switch (meta.rank) {
    #define BOOST_PP_LOCAL_MACRO(n) \
      case n: \
        if(meta.typeID=="CArray") result.extend(doRead<CArray<n>,n>(ifs,NPY_CDOUBLE)); \
        if(meta.typeID=="DArray") result.extend(doRead<DArray<n>,n>(ifs,NPY_DOUBLE));  \
        break;
    #define BOOST_PP_LOCAL_LIMITS (1, PYTHON_MAX_RANK)
    #include BOOST_PP_LOCAL_ITERATE()
  }
  ifs.close();
  return result;
}

void write(str filename, const numeric::array &array, double time)
{
  std::string f = extract<std::string>(filename);
  std::ofstream ofs(f.c_str(),std::ios_base::binary|std::ios_base::trunc);
  throw_file(ofs,f);

  PyArrayObject *a = reinterpret_cast<PyArrayObject *>(array.ptr());
  int rank=PyArray_NDIM(a);
  throw_rank(rank);
  switch (rank) {
    #define BOOST_PP_LOCAL_MACRO(n)                 \
      case n:                                       \
        if(PyArray_TYPE(a)==NPY_DOUBLE)            \
          doWrite<DArray<n>,double,n>(&ofs,a,time); \
        if(PyArray_TYPE(a)==NPY_CDOUBLE)           \
          doWrite<CArray<n>,dcomp,n>(&ofs,a,time);  \
        break;
    #define BOOST_PP_LOCAL_LIMITS (1, PYTHON_MAX_RANK)
    #include BOOST_PP_LOCAL_ITERATE()
  }
}


void pyDimensionsMismatchException(const trajectory::DimensionsMismatchException &)
{
  PyErr_SetString(PyExc_RuntimeError, "dimensions mismatch");
}

void pyRankMismatchException(const trajectory::RankMismatchException &)
{
  PyErr_SetString(PyExc_RuntimeError, "rank mismatch");
}

void export_io()
{
  import_array();
  numeric::array::set_module_and_type("numpy", "ndarray");
  def("read", read);
  def("write", write);

  register_exception_translator<trajectory::DimensionsMismatchException>(&pyDimensionsMismatchException);
  register_exception_translator<trajectory::RankMismatchException>      (&pyRankMismatchException);
  
  class_<trajectory::SerializationMetadata>("SerializationMetadata")
    .def_readonly("protocolVersion", &trajectory::SerializationMetadata::protocolVersion)
    .def_readonly("rank",            &trajectory::SerializationMetadata::rank)
    .def_readonly("trajectoryID",    &trajectory::SerializationMetadata::trajectoryID);
}

} // pythonext