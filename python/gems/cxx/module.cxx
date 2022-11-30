#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "pyKvlCalculator.h"
#include "pyKvlImage.h"
#include "pyKvlMesh.h"
#include "pyKvlOptimizer.h"
#include "pyKvlTransform.h"
#include "itkMultiThreader.h"

namespace py = pybind11;

void setGlobalDefaultNumberOfThreads(int maximumNumberOfThreads){
    itk::MultiThreader::SetGlobalDefaultNumberOfThreads( maximumNumberOfThreads );
}

PYBIND11_MODULE(gemsbindings, m) {
    py::class_<KvlImage>(m, "KvlImage")
            .def(py::init<const std::string &>())
            .def(py::init<const py::array_t<float> &>())
            .def(py::init<const std::string &, const std::string &>())
            .def_property_readonly("transform_matrix", &KvlImage::GetTransform)
            .def_property_readonly("non_cropped_image_size", &KvlImage::GetNonCroppedImageSize)
            .def_property_readonly("crop_slices", &KvlImage::GetCropSlices)
            .def("getImageBuffer", &KvlImage::GetImageBuffer)
            .def("write", &KvlImage::Write)
            .def("writeimage", &KvlImage::WriteImage)
            .def_static("smooth_image_buffer", &KvlImage::smoothImageBuffer)
            ;

    py::class_<KvlTransform>(m, "KvlTransform")
            .def(py::init<const py::array_t<double> &>())
            .def_property_readonly("as_numpy_array", &KvlTransform::AsNumpyArray)
            ;

    py::class_<KvlCostAndGradientCalculator>(m, "KvlCostAndGradientCalculator")
            .def(py::init<const std::string &,
                const std::vector<KvlImage> &,
                const std::string &,
                const KvlTransform &,
                const py::array_t<double> &,
                const py::array_t<double> &,
                const py::array_t<float> &,
                const py::array_t< int > &,
                const py::array_t<double> &>(),
                py::arg("typeName"),
                py::arg("images"),
                py::arg("boundaryCondition"),
                py::arg("transform")=KvlTransform(nullptr),
                py::arg("means")=py::array_t<double>(),
                py::arg("variances")=py::array_t<double>(),
                py::arg("mixtureWeights")=py::array_t<float>(),
                py::arg("numberOfGaussiansPerClass")=py::array_t<int>(),
                py::arg("targetPoints")=py::array_t<double>())
            .def(py::init<KvlMeshCollection &,
                const double &,
                const double &,
                const KvlTransform &>(),
                py::arg("meshCollection"),
                py::arg("K0"),
                py::arg("K1"),
                py::arg("transform"))
            .def("evaluate_mesh_position", &KvlCostAndGradientCalculator::EvaluateMeshPosition)
            // Aliases to help with profiling
            .def("evaluate_mesh_position_a", &KvlCostAndGradientCalculator::EvaluateMeshPosition)
            .def("evaluate_mesh_position_b", &KvlCostAndGradientCalculator::EvaluateMeshPosition)
            .def("evaluate_mesh_position_c", &KvlCostAndGradientCalculator::EvaluateMeshPosition)
            ;

    py::class_<KvlOptimizer>(m, "KvlOptimizer")
            .def(py::init<const std::string &, const KvlMesh&, const KvlCostAndGradientCalculator&, const std::map<std::string, double> &>())
            .def("step_optimizer", &KvlOptimizer::StepOptimizer)
            // Aliases to help with profiling
            .def("step_optimizer_warp", &KvlOptimizer::StepOptimizer)
            .def("step_optimizer_atlas", &KvlOptimizer::StepOptimizer)
            .def("step_optimizer_samseg", &KvlOptimizer::StepOptimizer)
            .def("update_calculator", &KvlOptimizer::UpdateCalculator)
            .def("update_mesh", &KvlOptimizer::UpdateMesh)
            ;

    py::class_<KvlMesh>(m, "KvlMesh")
            .def(py::init(), py::return_value_policy::take_ownership)
            .def_property_readonly("point_count", &KvlMesh::PointCount, py::return_value_policy::take_ownership)
            .def_property("points", &KvlMesh::GetPointSet, &KvlMesh::SetPointSet, py::return_value_policy::take_ownership)
            .def_property("alphas", &KvlMesh::GetAlphas, &KvlMesh::SetAlphas, py::return_value_policy::take_ownership)
            .def_property("can_moves", &KvlMesh::GetCanMoves, &KvlMesh::SetCanMoves, py::return_value_policy::take_ownership)
            .def( "fit_alphas", &KvlMesh::FitAlphas, py::return_value_policy::take_ownership )
            .def("draw_jacobian_determinant", &KvlMesh::DrawJacobianDeterminant, py::return_value_policy::take_ownership)
            .def("scale", &KvlMesh::Scale, py::return_value_policy::take_ownership)
            .def("rasterize_values", &KvlMesh::RasterizeValues, py::arg("shape"), py::arg("values"), py::return_value_policy::take_ownership)
            .def("rasterize", &KvlMesh::RasterizeMesh, py::arg("shape"), py::arg("classNumber") = -1, py::return_value_policy::take_ownership)
            .def("get_submesh", &KvlMesh::GetSubmesh, py::return_value_policy::take_ownership)
            // Aliases to help with profiling
            .def("rasterize_warp", &KvlMesh::RasterizeMesh, py::arg("shape"), py::arg("classNumber") = -1, py::return_value_policy::take_ownership)
            .def("rasterize_atlas", &KvlMesh::RasterizeMesh, py::arg("shape"), py::arg("classNumber") = -1, py::return_value_policy::take_ownership)
            .def("rasterize_1a", &KvlMesh::RasterizeMesh, py::arg("shape"), py::arg("classNumber") = -1, py::return_value_policy::take_ownership)
            .def("rasterize_1b", &KvlMesh::RasterizeMesh, py::arg("shape"), py::arg("classNumber") = -1, py::return_value_policy::take_ownership)
            .def("rasterize_2", &KvlMesh::RasterizeMesh, py::arg("shape"), py::arg("classNumber") = -1, py::return_value_policy::take_ownership)
            .def("rasterize_3", &KvlMesh::RasterizeMesh, py::arg("shape"), py::arg("classNumber") = -1, py::return_value_policy::take_ownership)
            ;

    py::class_<KvlMeshCollection>(m, "KvlMeshCollection")
            .def(py::init(), py::return_value_policy::take_ownership)
            .def_property_readonly("mesh_count", &KvlMeshCollection::MeshCount, py::return_value_policy::take_ownership)
            .def_property("k", &KvlMeshCollection::GetK, &KvlMeshCollection::SetK, py::return_value_policy::take_ownership)
            .def_property_readonly("reference_mesh", &KvlMeshCollection::GetReferenceMesh, py::return_value_policy::take_ownership)
            .def("get_mesh", &KvlMeshCollection::GetMesh, py::return_value_policy::take_ownership)
            .def_property("reference_position", &KvlMeshCollection::GetReferencePosition, &KvlMeshCollection::SetReferencePosition, py::return_value_policy::take_ownership)
            .def("set_positions", &KvlMeshCollection::SetPositions, py::return_value_policy::take_ownership)
            .def("construct", &KvlMeshCollection::Construct, py::return_value_policy::take_ownership)
            .def("read", &KvlMeshCollection::Read, py::return_value_policy::take_ownership)
            .def("transform", &KvlMeshCollection::Transform, py::return_value_policy::take_ownership)
            .def("smooth", &KvlMeshCollection::Smooth, py::return_value_policy::take_ownership)
            .def("write", &KvlMeshCollection::Write, py::return_value_policy::take_ownership)
            .def("generatefromsinglemesh", &KvlMeshCollection::GenerateFromSingleMesh, py::return_value_policy::take_ownership)
            ;
     m.def("setGlobalDefaultNumberOfThreads", &setGlobalDefaultNumberOfThreads, "Sets the maximum number of threads for ITK.");
}
