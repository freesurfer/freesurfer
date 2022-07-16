#pragma once

#include <algorithm>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;


template<class T>
using arrayf = py::array_t<T, py::array::forcecast | py::array::f_style>;

template<class T>
using arrayc = py::array_t<T, py::array::forcecast | py::array::c_style>;


enum MemoryOrder { C, Fortran };

std::vector<ssize_t> cstrides(std::vector<ssize_t> shape, size_t size);
std::vector<ssize_t> fstrides(std::vector<ssize_t> shape, size_t size);

std::string shapeString(std::vector<ssize_t> shape);

template<class T>
py::array_t<T> makeArray(std::vector<ssize_t> shape, std::vector<ssize_t> strides, const T* const data, bool free = true) {
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
  if (order == MemoryOrder::Fortran) {
    return makeArray(shape, fstrides(shape, sizeof(T)), data, free);
  } else {
    return makeArray(shape, cstrides(shape, sizeof(T)), data, free);
  }
}

template<class T>
py::array_t<T> copyArray(std::vector<ssize_t> shape, MemoryOrder order, const T* const data) {
  if (order == MemoryOrder::Fortran) {
    return py::array_t<T>(shape, fstrides(shape, sizeof(T)), data);
  } else {
    return py::array_t<T>(shape, cstrides(shape, sizeof(T)), data);
  }
}
