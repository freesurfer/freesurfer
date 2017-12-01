#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
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
            .def("getImageBuffer", &KvlImage::GetImageBuffer)
            ;

    py::class_<KvlCostAndGradientCalculator>(m, "KvlCostAndGradientCalculator")
            .def(py::init<const std::string &, const std::vector<KvlImage> &, const std::string &>())
            .def("evaluate_mesh_position", &KvlCostAndGradientCalculator::EvaluateMeshPosition)
            ;

    py::class_<KvlOptimizer>(m, "KvlOptimizer")
            .def(py::init<const std::string &, const KvlMesh&, const KvlCostAndGradientCalculator&, const std::map<std::string, double> &>())
            .def("step_optimizer", &KvlOptimizer::StepOptimizer)
            ;

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
            .def_property_readonly("reference_mesh", &KvlMeshCollection::GetReferenceMesh)
            .def("get_mesh", &KvlMeshCollection::GetMesh)
            .def("construct", &KvlMeshCollection::Construct)
            .def("read", &KvlMeshCollection::Read)
            .def("write", &KvlMeshCollection::Write)
            ;
}
