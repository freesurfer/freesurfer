#ifndef GEMS_PYKVLTRANSFORM_H
#define GEMS_PYKVLTRANSFORM_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

#include "kvlCroppedImageReader.h"

typedef kvl::CroppedImageReader::TransformType TransformType;
typedef TransformType::Pointer TransformPointer;


class KvlTransform {
public:
    // Python accessible
    KvlTransform(const py::array_t<double> &transformMatrix);
    py::array_t<double> AsNumpyArray() const;

    // C++ use only
    TransformPointer m_transform;

    KvlTransform(TransformPointer transform) : m_transform(transform) {};
    const TransformPointer GetTransform() const {
        return m_transform;
    }

};
#endif //GEMS_PYKVLTRANSFORM_H
