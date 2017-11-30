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
            .def("greet", &KvlImage::Greet)
            .def("getTransformMatrix", &KvlImage::GetTransformMatrix)
            ;

    m.def("useImage", &UseImage);

    py::class_<KvlMesh>(m, "KvlMesh")
            .def(py::init())
            ;

    py::class_<KvlMeshCollection>(m, "KvlMeshCollection")
            .def(py::init())
            .def("read", &KvlMeshCollection::Read)
            .def("write", &KvlMeshCollection::Write)
            .def("get_k", &KvlMeshCollection::GetK)
            .def("set_k", &KvlMeshCollection::SetK)
            .def("get_mesh", &KvlMeshCollection::GetMesh)
            .def("get_reference_mesh", &KvlMeshCollection::GetReferenceMesh)
            ;
}
