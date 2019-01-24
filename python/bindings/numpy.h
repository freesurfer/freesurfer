#ifndef NUMPY_H
#define NUMPY_H

#include <algorithm>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;


enum MemoryOrder {C, Fortran};

std::vector<ssize_t> cstrides(std::vector<ssize_t> shape, size_t size);
std::vector<ssize_t> fstrides(std::vector<ssize_t> shape, size_t size);

std::string shapeString(std::vector<ssize_t> shape);

template<class T>
py::array_t<T> makeArray(std::vector<ssize_t> shape, std::vector<ssize_t> strides, const T* const data, bool free = true) {
  // make python capsule handle
  py::capsule capsule;
  if (free) {
    capsule = py::capsule(data, [](void *d) { delete[] (T *)d; });
  } else {
    capsule = py::capsule(data);
  }
  return py::array_t<T>(shape, strides, data, capsule);
}

template<class T>
py::array_t<T> makeArray(std::vector<ssize_t> shape, MemoryOrder order, const T* const data, bool free = true) {
  // determine buffer array strides
  std::vector<ssize_t> strides;
  if (order == MemoryOrder::Fortran) {
    strides = fstrides(shape, sizeof(T));    
  } else {
    strides = cstrides(shape, sizeof(T));    
  }
  return makeArray(shape, strides, data, free);
}

#endif
