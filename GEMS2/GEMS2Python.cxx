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
            .def(py::init<const py::array_t<float> &>())
            .def(py::init<const std::string &, const std::string &>())
            .def_property_readonly("transform_matrix", &KvlImage::GetTransform)
            .def_property_readonly("non_cropped_image_size", &KvlImage::GetNonCroppedImageSize)
            .def_property_readonly("cropping_offset", &KvlImage::GetCroppingOffset)
            .def("getImageBuffer", &KvlImage::GetImageBuffer)
            .def("write", &KvlImage::Write)
            .def_static("smooth_image_buffer", &KvlImage::smoothImageBuffer)
            ;

    py::class_<KvlTransform>(m, "KvlTransform")
            .def(py::init<const py::array_t<double> &>())
            .def_property_readonly("as_numpy_array", &KvlTransform::AsNumpyArray)
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
            .def("scale", &KvlMesh::Scale)
            .def("rasterize", &KvlMesh::RasterizeMesh, py::arg("shape"), py::arg("classNumber") = -1)
            ;

    py::class_<KvlMeshCollection>(m, "KvlMeshCollection")
            .def(py::init())
            .def_property_readonly("mesh_count", &KvlMeshCollection::MeshCount)
            .def_property("k", &KvlMeshCollection::GetK, &KvlMeshCollection::SetK)
            .def_property_readonly("reference_mesh", &KvlMeshCollection::GetReferenceMesh)
            .def("get_mesh", &KvlMeshCollection::GetMesh)
            .def("construct", &KvlMeshCollection::Construct)
            .def("read", &KvlMeshCollection::Read)
            .def("transform", &KvlMeshCollection::Transform)
            .def("write", &KvlMeshCollection::Write)
            ;
}
