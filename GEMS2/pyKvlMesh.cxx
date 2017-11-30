#include <pybind11/pybind11.h>
#include "pyKvlMesh.h"
#include "itkMacro.h"

namespace py = pybind11;

KvlMeshCollection::KvlMeshCollection() {
    meshCollection = kvl::AtlasMeshCollection::New();
}

void KvlMeshCollection::Read(const std::string &meshCollectionFileName) {
    if (!meshCollection->Read(meshCollectionFileName.c_str())) {
        itkExceptionMacro("Couldn't read mesh collection from file " << meshCollectionFileName);
    }
    std::cout << "Read mesh collection: " << meshCollectionFileName << std::endl;
}

void KvlMeshCollection::Write(const std::string &meshCollectionFileName) {
    if (!meshCollection->Write(meshCollectionFileName.c_str())) {
        itkExceptionMacro("Couldn't write mesh collection to file " << meshCollectionFileName);
    }
    std::cout << "Wrote mesh collection: " << meshCollectionFileName << std::endl;
}

void KvlMeshCollection::SetK(double k) {
    meshCollection->SetK(k);
}

double KvlMeshCollection::GetK() const {
    return meshCollection->GetK();
}

KvlMesh::KvlMesh()
{

}
KvlMesh::KvlMesh(MeshPointer &aMesh) {
    mesh = aMesh;
}

KvlMesh *KvlMeshCollection::GetMesh(int meshNumber) {
    return NULL;
}

KvlMesh *KvlMeshCollection::GetReferenceMesh() {
    return NULL;
}
