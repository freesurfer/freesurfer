#ifndef GEMS_PYKVLMESH_H_H
#define GEMS_PYKVLMESH_H_H

#include "itkObject.h"
#include "kvlAtlasMeshCollection.h"

typedef kvl::AtlasMeshCollection::Pointer MeshCollectionPointer;
typedef kvl::AtlasMesh::ConstPointer MeshPointer;

class KvlMesh {
    MeshPointer mesh;
public:
    // Python accessible
    KvlMesh();

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

#endif //GEMS_PYKVLMESH_H_H
