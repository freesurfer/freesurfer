#ifndef GEMS_PYKVLTRANSFORM_H
#define GEMS_PYKVLTRANSFORM_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

#include "kvlCroppedImageReader.h"

typedef kvl::CroppedImageReader::TransformType TransformType;
typedef TransformType::Pointer TransformPointer;

py::array_t<double> TransformToNumpy(TransformPointer transform);
TransformPointer NumpyToTransform(py::array_t<double> transform);

class KvlTransform {
    TransformPointer m_transform;
public:
    // Python accessible
    KvlTransform(const py::array_t<double> &transformMatrix);
    py::array_t<double> GetTransformMatrix() const;

    // C++ use only
    const TransformPointer GetTransform() {
        return m_transform;
    }
};
#endif //GEMS_PYKVLTRANSFORM_H
