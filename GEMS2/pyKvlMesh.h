#ifndef GEMS_PYKVLMESH_H_H
#define GEMS_PYKVLMESH_H_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <vector>
#include "itkObject.h"
#include "kvlAtlasMeshCollection.h"


namespace py = pybind11;

typedef kvl::AtlasMeshCollection::Pointer MeshCollectionPointer;
typedef kvl::AtlasMesh::ConstPointer MeshPointer;
typedef kvl::AtlasMesh::PointsContainer::ConstPointer PointSetConstPointer;
typedef kvl::AtlasMesh::PointDataContainer::ConstPointer PointDataConstPointer;
typedef std::vector<unsigned int> SHAPE_3D;

class KvlMesh {

public:
    // Python accessible
    KvlMesh();
    int PointCount() const;
    py::array_t<double> GetPointSet() const;
    void SetPointSet(py::array_t<double>);
    py::array_t<double> GetAlphas() const;
    void SetAlphas(py::array_t<double>);

    // C++ Only
    KvlMesh(MeshPointer& aMesh);
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
    void Construct(const SHAPE_3D &meshSize, const SHAPE_3D &domainSize,
                   double initialStiffness,
                   unsigned int numberOfClasses, unsigned int numberOfMeshes);

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
py::array_t<double> AlphasToNumpy(PointDataConstPointer alphas);

#endif //GEMS_PYKVLMESH_H_H
