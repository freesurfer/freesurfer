#ifndef GEMS_PYKVLMESH_H_H
#define GEMS_PYKVLMESH_H_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "itkObject.h"
#include "kvlAtlasMeshCollection.h"

namespace py = pybind11;

typedef kvl::AtlasMeshCollection::Pointer MeshCollectionPointer;
typedef kvl::AtlasMesh::ConstPointer MeshPointer;
typedef kvl::AtlasMesh::PointsContainer* PointSetPointer;

class KvlMesh {
    MeshPointer mesh;
public:
    // Python accessible
    KvlMesh();
    int PointCount() const;

    // C++ Only
    KvlMesh(MeshPointer& aMesh);
};

class KvlMeshCollection {
    MeshCollectionPointer  meshCollection;
public:
    // Python accessible
    KvlMeshCollection();
    void Read(const std::string &meshCollectionFileName);
    void Write(const std::string &meshCollectionFileName);
    double GetK() const;
    void SetK(double k);
    unsigned int MeshCount() const;
    KvlMesh* GetMesh(int meshNumber);
    KvlMesh* GetReferenceMesh();

    // C++ use only
    const char *GetNameOfClass() const {
        return "KvlMeshCollection";
    }
    MeshCollectionPointer GetMeshCollection() {
        return meshCollection;
    }
};

py::array_t<double> PointSetToNumpy(PointSetPointer points);

#endif //GEMS_PYKVLMESH_H_H
