#ifndef GEMS_PYKVLNUMPY_H
#define GEMS_PYKVLNUMPY_H

#include <algorithm>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;


template<class T>
py::array_t<T> createNumpyArray(std::vector<size_t> shape, std::vector<size_t> strides, T *data){
    py::capsule free_when_done(data, [](void *d) { delete[] (T *)d; });
    auto result = py::array_t<T>(shape, strides, data, free_when_done);
    return result;
}

template<class T>
py::array_t<T> createNumpyArrayCStyle(std::vector<size_t> shape, const T* const data){
    size_t size = sizeof(T);
    std::vector<size_t> strides = {1};
    for(size_t d = shape.size()-1; d > 0; d--){
        strides.push_back(strides.back()*shape[d]);
    }
    std::reverse(strides.begin(), strides.end());
    for(auto &stride: strides){
        stride *= size;
    }
    return createNumpyArray(shape, strides, data);
}

template<class T>
py::array_t<T> createNumpyArrayFStyle(std::vector<size_t> shape, const T* const data){
    size_t size = sizeof(T);
    std::vector<size_t> strides = {1};
    for(size_t d = 0; d < shape.size()-1; d++){
        strides.push_back(strides.back()*shape[d]);
    }
    for(auto &stride: strides){
        stride *= size;
    }
    return createNumpyArray(shape, strides, data);
}

#endif //GEMS_PYKVLNUMPY_H
