#include <pybind11/pybind11.h>
#include "pyKvlCalculator.h"
#include "pyKvlImage.h"
#include "pyKvlMesh.h"
#include "pyKvlOptimizer.h"
#include "pyKvlTransform.h"

namespace py = pybind11;

PYBIND11_MODULE(GEMS2Python, m) {
    py::class_<KvlImage>(m, "KvlImage")
            .def(py::init<const std::string &>())
            .def("greet", &KvlImage::greet);

    m.def("useImage", &useImage);

    py::class_<KvlMeshCollection>(m, "KvlMeshCollection")
            .def(py::init<const std::string &>());

}
