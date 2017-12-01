#include "kvlAtlasMesh.h"
#include <pybind11/pybind11.h>
#include "itkMacro.h"
#include "pyKvlMesh.h"
#include "pyKvlNumpy.h"

KvlMesh::KvlMesh() {

}

KvlMesh::KvlMesh(MeshPointer &aMesh) {
    mesh = aMesh;
}

int KvlMesh::PointCount() const {
    return mesh->GetPoints()->Size();
}

py::array_t<double> KvlMesh::GetPointSet() const {
    return PointSetToNumpy(mesh->GetPoints());
}

py::array_t<double> KvlMesh::GetAlphas() const {
    return AlphasToNumpy(mesh->GetPointData());
}

void KvlMesh::SetPointSet(const py::array_t<double> &source) {
    PointSetPointer points = const_cast<PointSetPointer>(mesh->GetPoints());
    CopyNumpyToPointSet(points, source);
}

void KvlMesh::SetAlphas(const py::array_t<double> &source) {
    PointDataPointer alphas = const_cast<PointDataPointer>(mesh->GetPointData());
    CopyNumpyToPointDataSet(alphas, source);
}

KvlMeshCollection::KvlMeshCollection() {
    meshCollection = kvl::AtlasMeshCollection::New();
}

void KvlMeshCollection::Construct(const SHAPE_3D &meshSize, const SHAPE_3D &domainSize,
                                  double initialStiffness,
                                  unsigned int numberOfClasses, unsigned int numberOfMeshes) {
    if (meshSize.size() != 3) {
        itkExceptionMacro("meshSize must have 3 values not " << meshSize.size());
    }
    if (domainSize.size() != 3) {
        itkExceptionMacro("domainSize must have 3 values not " << domainSize.size());
    }
    unsigned int meshShape[3];
    unsigned int domainShape[3];
    for (int index = 0; index < 3; index++) {
        meshShape[index] = meshSize[index];
        domainShape[index] = domainSize[index];
    }
    meshCollection->Construct(meshShape, domainShape, initialStiffness,
                              numberOfClasses, numberOfMeshes);
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


KvlMesh *KvlMeshCollection::GetMesh(int meshNumber) {
    if (meshNumber < -1) {
        itkExceptionMacro("meshNumber " << meshNumber << " too negative");
    } else if (meshNumber == -1) {
        return this->GetReferenceMesh();
    } else {
        unsigned int meshCount = this->MeshCount();
        unsigned int meshIndex = (unsigned int) meshNumber;
        if (meshIndex >= meshCount) {
            itkExceptionMacro("meshNumber " << meshNumber << " exceeds mesh count of " << meshCount);
        } else {
            MeshPointer mesh = meshCollection->GetMesh(meshIndex);
            return new KvlMesh(mesh);
        }
    }
}

KvlMesh *KvlMeshCollection::GetReferenceMesh() {
    MeshPointer mesh = meshCollection->GetReferenceMesh();
    return new KvlMesh(mesh);
}

unsigned int KvlMeshCollection::MeshCount() const {
    return meshCollection->GetNumberOfMeshes();
}

#define XYZ_DIMENSIONS 3

py::array_t<double> PointSetToNumpy(PointSetConstPointer points) {
    const unsigned int numberOfNodes = points->Size();
    auto *data = new double[numberOfNodes * XYZ_DIMENSIONS];
    auto dataIterator = data;
    for (auto pointsIterator = points->Begin(); pointsIterator != points->End(); ++pointsIterator) {
        for (int xyzAxisSelector = 0; xyzAxisSelector < XYZ_DIMENSIONS; xyzAxisSelector++) {
            *dataIterator++ = pointsIterator.Value()[xyzAxisSelector];
        }
    }
    return createNumpyArrayCStyle({numberOfNodes, XYZ_DIMENSIONS}, data);
}

void CopyNumpyToPointSet(PointSetPointer points, const py::array_t<double> &source) {
    if (source.ndim() != 2) {
        itkGenericExceptionMacro("point source shape must have two dimensions");
    }
    const unsigned int currentNumberOfNodes = points->Size();
    auto sourceShape = source.shape();
    const unsigned int sourceNumberOfNodes = *sourceShape++;
    const unsigned int sourceXYZDimensions = *sourceShape++;
    if (sourceXYZDimensions != 3) {
        itkGenericExceptionMacro("source points must have 3 coordinates not " << sourceXYZDimensions);
    }
    if (sourceNumberOfNodes != currentNumberOfNodes) {
        itkGenericExceptionMacro(
                "source point count of "
                        << sourceNumberOfNodes
                        << " not equal to mesh point count of "
                        << currentNumberOfNodes
        );
    }
    unsigned int pointIndex = 0;
    for (auto pointsIterator = points->Begin(); pointsIterator != points->End(); ++pointsIterator, ++pointIndex) {
        for (int xyzAxisSelector = 0; xyzAxisSelector < XYZ_DIMENSIONS; xyzAxisSelector++) {
            pointsIterator.Value()[xyzAxisSelector] = *source.data(pointIndex, xyzAxisSelector);
        }
    }
}

py::array_t<double> AlphasToNumpy(PointDataConstPointer alphas) {
    const unsigned int numberOfNodes = alphas->Size();
    const unsigned int numberOfLabels = alphas->Begin().Value().m_Alphas.Size();
    auto *data = new double[numberOfNodes * numberOfLabels];
    auto dataIterator = data;
    for (auto alphasIterator = alphas->Begin(); alphasIterator != alphas->End(); ++alphasIterator) {
        for (int label = 0; label < numberOfLabels; label++) {
            *dataIterator++ = alphasIterator.Value().m_Alphas[label];
        }
    }
    return createNumpyArrayCStyle({numberOfNodes, numberOfLabels}, data);
}

void CopyNumpyToPointDataSet(PointDataPointer destinationAlphas, const py::array_t<double> &source)
{
    if (source.ndim() != 2) {
        itkGenericExceptionMacro("data point source shape must have two dimensions");
    }
    const unsigned int currentNumberOfNodes = destinationAlphas->Size();
    auto sourceShape = source.shape();
    const unsigned int sourceNumberOfNodes = *sourceShape++;
    const unsigned int numberOfLabels = *sourceShape++;
    if (sourceNumberOfNodes != currentNumberOfNodes) {
        itkGenericExceptionMacro(
                "source data point count of "
                        << sourceNumberOfNodes
                        << " not equal to mesh point count of "
                        << currentNumberOfNodes
        );
    }
    if (numberOfLabels <= 0) {
        itkGenericExceptionMacro("source data have positive number of labels not " << numberOfLabels);
    }
    unsigned int pointIndex = 0;
    for (auto alphasIterator = destinationAlphas->Begin(); alphasIterator != destinationAlphas->End(); ++alphasIterator, ++pointIndex) {
        kvl::AtlasAlphasType  alphas( numberOfLabels );
        for (int label = 0; label < numberOfLabels; label++) {
            alphas[label] = *source.data(pointIndex, label);
        }
        alphasIterator.Value().m_Alphas = alphas;
    }
}
