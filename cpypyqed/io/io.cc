// -*- C++ -*-
#include "Evolved.h"
#include "PythonExtension.h"
#include "Trajectory.tcc"
#include "Types.h"

#include <blitz/tinyvec2.h>
#include <boost/function.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/preprocessor/iteration/local.hpp>
#include <numpy/ndarrayobject.h>
#include <fstream>
#include <string>
#include <iostream>

using namespace trajectory;
using namespace boost::python;

namespace pythonext {

namespace {

template<int RANK>
object doRead(std::ifstream &ifs)
{
  list states;
  list times;
  list result;
  CArray<RANK> a(1);
  a=0;
  AdaptiveIO<CArray<RANK>> traj(evolved::makeIO(a));
  while ( (ifs.peek(), !ifs.eof()) ) {
    trajectory::readViaSStream(traj,ifs);
    npy_intp dims[RANK];
    for (int i=0;i<RANK;i++) dims[i]=a.extent(i);
    PyObject * pyObj = PyArray_SimpleNewFromData(RANK, dims, NPY_CDOUBLE,a.dataFirst());
    handle<> h( pyObj );
    numeric::array arr( h );
    states.append(arr.copy());
    times.append(traj.getEvolvedIO()->getTime());
  }
  return make_tuple(states,times);
}

template<int RANK>
void doWrite(std::ofstream *ofs, const PyArrayObject *a)
{
  npy_intp *dims=PyArray_DIMS(a);
  blitz::TinyVector<int,RANK> shape;
  for (int i=0; i<RANK; i++) shape[i]=dims[i];
  CArray<RANK>  blitz_a = CArray<RANK>(static_cast<dcomp *>(PyArray_DATA(a)), shape, blitz::duplicateData);
  AdaptiveIO<CArray<RANK>> traj(evolved::makeIO(blitz_a));
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

}

object read(str filename) 
{
  std::string f = extract<std::string>(filename);

  std::ifstream ifs(f.c_str(),std::ios_base::binary);
  throw_file(ifs,f);

  trajectory::SerializationMetadata meta = trajectory::readMeta(ifs);
  
//   if (!(meta.trajectoryType=="Master"||meta.trajectoryType=="MCWF_Trajectory"||meta.trajectoryType=="Dummy") ){
//     PyErr_SetString(PyExc_NotImplementedError, (meta.trajectoryType+std::string(": not implemented")).c_str());
//     throw_error_already_set();
//   }
  throw_rank(meta.rank);
  
  list result;
  result.append(meta);
  switch (meta.rank) {
    #define BOOST_PP_LOCAL_MACRO(n) \
      case n: result.extend(doRead<n>(ifs)); break;
    #define BOOST_PP_LOCAL_LIMITS (1, PYTHON_MAX_RANK)
    #include BOOST_PP_LOCAL_ITERATE()
  }
  ifs.close();
  return result;
}

void write(const numeric::array &array, str filename)
{
  std::string f = extract<std::string>(filename);
  std::ofstream ofs(f.c_str(),std::ios_base::binary|std::ios_base::trunc);
  throw_file(ofs,f);

  const PyArrayObject *a = reinterpret_cast<const PyArrayObject *>(array.ptr());
  int rank=PyArray_NDIM(a);
  throw_rank(rank);

  switch (rank) {
    #define BOOST_PP_LOCAL_MACRO(n)               \
      case n: doWrite<n>(&ofs,a); break;
    #define BOOST_PP_LOCAL_LIMITS (1, PYTHON_MAX_RANK)
    #include BOOST_PP_LOCAL_ITERATE()
  }
}

void export_io()
{
  import_array();
  numeric::array::set_module_and_type("numpy", "ndarray");
  def("read", read);
  def("write", write);
  
  class_<trajectory::SerializationMetadata>("SerializationMetadata")
    .def_readonly("protocolVersion", &trajectory::SerializationMetadata::protocolVersion)
    .def_readonly("rank",            &trajectory::SerializationMetadata::rank)
    .def_readonly("trajectoryID",    &trajectory::SerializationMetadata::trajectoryID);
}

} // pythonext