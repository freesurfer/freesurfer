#include <sstream>

#include "bindings_numpy.h"


/*
  Computes C-contiguous strides for a buffer of shape `shape` and element size `size`.
*/
std::vector<ssize_t> cstrides(std::vector<ssize_t> shape, size_t size) {
  std::vector<ssize_t> strides = {1};
  for(unsigned int d = shape.size()-1; d > 0; d--) {
    strides.push_back(strides.back() * shape[d]);
  }
  std::reverse(strides.begin(), strides.end());
  for (ssize_t &stride : strides) {
    stride *= size;
  }
  return strides;
}


/*
  Computes Fortran-contiguous strides for a buffer of shape `shape` and element size `size`.
*/
std::vector<ssize_t> fstrides(std::vector<ssize_t> shape, size_t size) {
  std::vector<ssize_t> strides = {1};
  for (unsigned int d = 0; d < shape.size()-1; d++) {
    strides.push_back(strides.back() * shape[d]);
  }
  for (ssize_t &stride : strides) {
    stride *= size;
  }
  return strides;
}


/*
  Converts a shape vector into a readable string, such as (256, 256, 256, 1).
*/
std::string shapeString(std::vector<ssize_t> shape){
  std::stringstream ss;
  ss << "(";
  for (unsigned int i = 0; i < shape.size(); i++) {
    if (i != 0) ss << ", ";
    ss << shape[i];
  }
  ss << ")";
  return ss.str();
}
