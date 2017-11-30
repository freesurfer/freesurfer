#ifndef GEMS_PYKVLNUMPY_H
#define GEMS_PYKVLNUMPY_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

template<class T>
py::array_t<T> createNumpyArray(std::vector<size_t> shape, T *data){
    size_t total_size = 1;
    for( auto &i : shape) {
        total_size *= i;
    }
    py::capsule free_when_done(data, [](void *d) {
        std::cout << "destructor was called on numpy data";
        delete[] (T *)d;
    });
    auto result = py::array_t<T>(total_size, data, free_when_done);
    result.resize(shape);
    return result;
}

template<class T>
py::array_t<T> createNumpyArray(std::vector<size_t> shape){
    size_t total_size = 1;
    for( auto &i : shape) {
        total_size *= i;
    }
    auto *data = new T[total_size];
    return createNumpyArray(shape, data);
}

#endif //GEMS_PYKVLNUMPY_H
