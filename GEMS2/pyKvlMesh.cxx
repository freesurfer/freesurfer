#include "kvlAtlasMesh.h"
#include <pybind11/pybind11.h>
#include "itkMacro.h"
#include "pyKvlMesh.h"
#include "pyKvlTransform.h"
#include "pyKvlNumpy.h"

#define XYZ_DIMENSIONS 3

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

void TransformMeshPoints(MeshPointer mesh, const double *scaleFactor) {
    PointSetPointer points = const_cast<PointSetPointer>(mesh->GetPoints());
    for (auto pointsIterator = points->Begin(); pointsIterator != points->End(); ++pointsIterator) {
        for (int xyzIndex = 0; xyzIndex < XYZ_DIMENSIONS; xyzIndex++) {
            pointsIterator.Value()[xyzIndex] *= scaleFactor[xyzIndex];
        }
    }
}

void TransformCellData(MeshPointer mesh, const double *scaling) {
    // Also the reference position of the mesh has changed. Note, however, that we don't
    // actually have access to the original reference position, only some sufficient
    // statistics calculated from it. In particular, we have only the three first columns
    // of the matrix Z = inv( X ) = inv( [ p0 p1 p2 p3; 1 1 1 1 ] ). The transformed
    // position Xtrans is given by
    //
    //     Xtrans = diag( scaling[0] scaling[ 1 ] scaling[ 2 ] 1 ) * X + [ translation[ 0 ] translation[ 1 ] translation[ 2 ] 1 ]'
    //
    // but since we'll also end up calculating the upper 3x3 matrix of Ytrans * Ztrans
    // (with Y the equivalent of X but in the deformed mesh)
    // to evaluate the deformation penatly, we can safely drop the translation part.
    // In short, we will calculate the first 3 columns of the matrix
    //
    //    Ztrans = inv( diag( scaling[0] scaling[ 1 ] scaling[ 2 ] 1 ) * X )
    //
    //           = Z * diag( 1/scaling[0] 1/scaling[1] 1/scaling[2] 1 )
    //
    // which is given by multiplying each column i of Z with a factor 1/scaling[i]
    //
    kvl::AtlasMesh::CellDataContainer* cellData = const_cast<kvl::AtlasMesh::CellDataContainer*>(mesh->GetCellData());
    for ( kvl::AtlasMesh::CellDataContainer::Iterator  it = cellData->Begin();
          it != mesh->GetCellData()->End(); ++it )
    {
        it.Value().m_ReferenceVolumeTimesK *= ( scaling[ 0 ] * scaling[ 1 ] * scaling[ 2 ] );

        it.Value().m_Z11 /= scaling[ 0 ];
        it.Value().m_Z21 /= scaling[ 0 ];
        it.Value().m_Z31 /= scaling[ 0 ];
        it.Value().m_Z41 /= scaling[ 0 ];

        it.Value().m_Z12 /= scaling[ 1 ];
        it.Value().m_Z22 /= scaling[ 1 ];
        it.Value().m_Z32 /= scaling[ 1 ];
        it.Value().m_Z42 /= scaling[ 1 ];

        it.Value().m_Z13 /= scaling[ 2 ];
        it.Value().m_Z23 /= scaling[ 2 ];
        it.Value().m_Z33 /= scaling[ 2 ];
        it.Value().m_Z43 /= scaling[ 2 ];
    }
}

void KvlMesh::Scale(const SCALE_3D &scaling) {
    auto scaleShape = scaling.size();
    double scaleFactor[XYZ_DIMENSIONS];
    if(scaleShape == 1) {
        scaleFactor[0] = scaleFactor[1] = scaleFactor[2] = scaling[0];
    } else if (scaleShape == 3) {
        for(int xyzIndex=0; xyzIndex < XYZ_DIMENSIONS; xyzIndex++) {
            scaleFactor[xyzIndex] = scaling[xyzIndex];
        }
    } else {
        itkExceptionMacro("mesh scaling dimensions must be either 1 or 3 not " << scaleShape);
    }
    TransformMeshPoints(mesh, scaleFactor);
    TransformCellData(mesh, scaleFactor);
}

void KvlMesh::Transform(const py::array_t<double> &transformMatrix) {
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

void KvlMeshCollection::Transform(const py::array_t<double> &transformData) {
    TransformPointer transform = NumpyTo
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
