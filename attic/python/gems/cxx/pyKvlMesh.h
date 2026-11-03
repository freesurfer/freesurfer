#ifndef GEMS_PYKVLMESH_H_H
#define GEMS_PYKVLMESH_H_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <vector>
#include "itkObject.h"
#include "kvlAtlasMeshCollection.h"
#include "kvlAtlasMeshSmoother.h"
#include "pyKvlTransform.h"

namespace py = pybind11;

typedef kvl::AtlasMeshCollection::Pointer MeshCollectionPointer;
typedef kvl::AtlasMesh::ConstPointer MeshPointer;
typedef kvl::AtlasMesh::PointsContainer* PointSetPointer;
typedef kvl::AtlasMesh::PointsContainer::ConstPointer PointSetConstPointer;
typedef kvl::AtlasMesh::PointDataContainer* PointDataPointer;
typedef kvl::AtlasMesh::PointDataContainer::ConstPointer PointDataConstPointer;
typedef std::vector<unsigned int> SHAPE_3D;
typedef std::vector<double> SCALE_3D;

class KvlMesh {

public:
    // Python accessible
    KvlMesh();
    int PointCount() const;
    py::array_t<double> GetPointSet() const;
    void SetPointSet(const py::array_t<double> &source);
    py::array_t<double> GetAlphas() const;
    void SetAlphas(const py::array_t<double> &source);
    py::array_t<bool> GetCanMoves() const;
    void SetCanMoves(const py::array_t<bool> &source);
    void Scale(const SCALE_3D &scaling);
    py::array_t<uint16_t> RasterizeMesh(std::vector<size_t> size, int classNumber=-1);
    py::array RasterizeValues(std::vector<size_t> size, py::array_t<double, py::array::c_style | py::array::forcecast> values);
    py::array_t<double> FitAlphas( const py::array_t< uint16_t, py::array::f_style | py::array::forcecast >& probabilityImageBuffer, int EMIterations=10 ) const;
    py::array_t<double> DrawJacobianDeterminant(std::vector<size_t> size);
    KvlMesh* GetSubmesh( py::array_t<bool>& mask );

    // C++ Only
    KvlMesh(MeshPointer& aMesh);
    const char *GetNameOfClass() const {
        return "KvlMesh";
    }
    MeshPointer mesh;
};

class KvlMeshCollection {
    MeshCollectionPointer  meshCollection;
public:
    // Python accessible
    void Read(const std::string &meshCollectionFileName);
    void Write(const std::string &meshCollectionFileName);
    double GetK() const;
    void SetK(double k);
    unsigned int MeshCount() const;
    KvlMesh* GetMesh(int meshNumber);
    KvlMesh* GetReferenceMesh();
    py::array_t<double>  GetReferencePosition() const;
    void SetReferencePosition(const py::array_t<double> &source);
    void SetPositions(const py::array_t<double> &reference, const std::vector<py::array_t<double>> &positions);
    void Construct(const SHAPE_3D &meshSize, const SHAPE_3D &domainSize,
                   double initialStiffness,
                   unsigned int numberOfClasses, unsigned int numberOfMeshes);
    void Transform(const KvlTransform &transform);
    void Smooth(double sigma);
    void GenerateFromSingleMesh(const KvlMesh &singleMesh, unsigned int numberOfMeshes, double K);

    // C++ use only
    KvlMeshCollection();
    const char *GetNameOfClass() const {
        return "KvlMeshCollection";
    }
    MeshCollectionPointer GetMeshCollection() {
        return meshCollection;
    }
};

py::array_t<double> PointSetToNumpy(PointSetConstPointer points);
void CopyNumpyToPointSet(PointSetPointer points, const py::array_t<double> &source);
py::array_t<double> AlphasToNumpy(PointDataConstPointer alphas);
py::array_t<bool> CanMovesToNumpy(PointDataConstPointer pointData );
void CopyNumpyToPointDataSet(PointDataPointer alphas, const py::array_t<double> &source);
void CopyNumpyToCanMoves(PointDataPointer destinationPointData, const py::array_t<bool> &source);
void CreatePointSetFromNumpy(PointSetPointer targetPoints, const py::array_t<double> &source);

#endif //GEMS_PYKVLMESH_H_H
