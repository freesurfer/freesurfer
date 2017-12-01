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
            .def("getTransformMatrix", &KvlImage::GetTransformMatrix)
            .def("getImageBuffer", &KvlImage::GetImageBuffer);

    py::class_<KvlMesh>(m, "KvlMesh")
            .def(py::init())
            .def_property_readonly("point_count", &KvlMesh::PointCount)
            .def_property("points", &KvlMesh::GetPointSet, &KvlMesh::SetPointSet)
            .def_property("alphas", &KvlMesh::GetAlphas, &KvlMesh::SetAlphas)
            ;

    py::class_<KvlMeshCollection>(m, "KvlMeshCollection")
            .def(py::init())
            .def_property_readonly("mesh_count", &KvlMeshCollection::MeshCount)
            .def_property("k", &KvlMeshCollection::GetK, &KvlMeshCollection::SetK)
            .def("get_reference_mesh", &KvlMeshCollection::GetReferenceMesh)
            .def("get_mesh", &KvlMeshCollection::GetMesh)
            .def("construct", &KvlMeshCollection::Construct)
            .def("read", &KvlMeshCollection::Read)
            .def("write", &KvlMeshCollection::Write)
            ;
}
