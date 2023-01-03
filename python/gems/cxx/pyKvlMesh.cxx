#include "kvlAtlasMesh.h"
#include <pybind11/pybind11.h>
#include <itkImageRegionConstIterator.h>
#include "itkMacro.h"
#include "pyKvlMesh.h"
#include "pyKvlTransform.h"
#include "pyKvlNumpy.h"
#include "vnl/vnl_det.h"
#include "itkCellInterface.h"
#include "kvlAtlasMeshAlphaDrawer.h"
#include "kvlAtlasMeshMultiAlphaDrawer.h"
#include "kvlAtlasMeshValueDrawer.h"
#include "kvlAtlasMeshProbabilityImageStatisticsCollector.h"
#include "kvlAtlasMeshJacobianDeterminantDrawer.h"
#include "itkImageRegionIterator.h"


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

py::array_t<bool> KvlMesh::GetCanMoves() const {
    return CanMovesToNumpy(mesh->GetPointData());
}

void KvlMesh::SetPointSet(const py::array_t<double> &source) {
    PointSetPointer points = const_cast<PointSetPointer>(mesh->GetPoints());
    CopyNumpyToPointSet(points, source);
}

void KvlMesh::SetAlphas(const py::array_t<double> &source) {
    PointDataPointer alphas = const_cast<PointDataPointer>(mesh->GetPointData());
    CopyNumpyToPointDataSet(alphas, source);
}

void KvlMesh::SetCanMoves(const py::array_t<bool> &source) {
    PointDataPointer pointData = const_cast<PointDataPointer>(mesh->GetPointData());
    CopyNumpyToCanMoves(pointData, source);
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

void ApplyTransformToMeshes(MeshCollectionPointer meshCollection, const TransformPointer transform)
{
    meshCollection->Transform( -1, transform );
    for ( unsigned int i = 0; i < meshCollection->GetNumberOfMeshes(); i++ )
    {
        meshCollection->Transform( i, transform );
    }
}

void ReverseTetrahedronSidedness(MeshCollectionPointer meshCollection) {
    //std::cout << "Careful here: the applied transformation will turn positive tetrahedra into negative ones." << std::endl;
    //std::cout << transform->GetMatrix().GetVnlMatrix() << std::endl;
    //std::cout << " determinant: " << determinant << std::endl;
    //std::cout << "Starting to swap the point assignments of each tetrahedron..." << std::endl;
    kvl::AtlasMesh::CellsContainer* cells = meshCollection-> GetCells();

    for ( kvl::AtlasMesh::CellsContainer::Iterator  cellIt = cells->Begin();
          cellIt != meshCollection->GetCells()->End(); ++cellIt )
    {
        kvl::AtlasMesh::CellType*  cell = cellIt.Value();

        if( cell->GetType() != kvl::AtlasMesh::CellType::TETRAHEDRON_CELL )
        {
            continue;
        }

        // Swap points assigned to first two vertices. This will readily turn negative tetrahedra
        //into positives ones.
        kvl::AtlasMesh::CellType::PointIdIterator  pit = cell->PointIdsBegin();
        const kvl::AtlasMesh::PointIdentifier  p0Id = *pit;
        ++pit;
        const kvl::AtlasMesh::PointIdentifier  p1Id = *pit;

        pit = cell->PointIdsBegin();
        *pit = p1Id;
        ++pit;
        *pit = p0Id;
    } // End loop over all tetrahedra

}

void KvlMeshCollection::Transform(const KvlTransform &kvlTransform) {
    auto transform = kvlTransform.GetTransform();
    ApplyTransformToMeshes(meshCollection, transform);
    const float  determinant = vnl_det( transform->GetMatrix().GetVnlMatrix() );
    if (determinant < 0) {
        ReverseTetrahedronSidedness(meshCollection);
    }
}


void KvlMeshCollection::Write(const std::string &meshCollectionFileName) {
    if (!meshCollection->Write(meshCollectionFileName.c_str())) {
        itkExceptionMacro("Couldn't write mesh collection to file " << meshCollectionFileName);
    }
    std::cout << "Wrote mesh collection: " << meshCollectionFileName << std::endl;
}

void KvlMeshCollection::GenerateFromSingleMesh(const KvlMesh &singleMesh, unsigned int numberOfMeshes, double K)
{
  kvl::AtlasMesh::ConstPointer constMesh = singleMesh.mesh;
  kvl::AtlasMesh::Pointer mutableMesh = const_cast< kvl::AtlasMesh* >( constMesh.GetPointer() );
  meshCollection->GenerateFromSingleMesh(mutableMesh, numberOfMeshes, K);
}

void KvlMeshCollection::SetK(double k) {
    meshCollection->SetK(k);
}

double KvlMeshCollection::GetK() const {
    return meshCollection->GetK();
}

py::array_t<double> KvlMeshCollection::GetReferencePosition() const {
  //KvlMesh mesh = KvlMeshCollection::GetMesh(-1);
  return PointSetToNumpy(meshCollection->GetReferencePosition());
}

void KvlMeshCollection::SetReferencePosition(const py::array_t<double> &source) {
  PointSetPointer points = const_cast<PointSetPointer>(meshCollection->GetReferencePosition());
  CopyNumpyToPointSet(points, source);
  meshCollection->SetReferencePosition(points);
}

void KvlMeshCollection::SetPositions(const py::array_t<double> &reference, const std::vector<py::array_t<double>> &positions) {

  // set the reference position
  SetReferencePosition(reference);
  std::cout << "positions.size(): " << positions.size() << std::endl;
  // set the additional positions
  std::vector<kvl::AtlasMeshCollection::PointsContainerType::Pointer>  pointsVector;
  for (auto &pos : positions) {
    // Create a new Points container
    typedef kvl::AtlasMesh::PointsContainer  PointsContainerType;
    PointsContainerType::Pointer  target = PointsContainerType::New();
    PointsContainerType::ConstPointer  sourcePosition = meshCollection->GetReferencePosition();
    PointsContainerType::ConstIterator  sourceIt = sourcePosition->Begin();
    while ( sourceIt != sourcePosition->End() )
      {
      // insert source coords for this point into target
      target->InsertElement( sourceIt.Index(), sourceIt.Value() );
      ++sourceIt;
      }

    CopyNumpyToPointSet( target.GetPointer(), pos);
    pointsVector.push_back( target );
  }
  meshCollection->SetPositions(pointsVector);
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

void KvlMeshCollection::Smooth(double sigma) {
    kvl::AtlasMeshCollection::Pointer atlasMeshCollectionPtr = this->GetMeshCollection();
    kvl::AtlasMeshSmoother::Pointer  smoother = kvl::AtlasMeshSmoother::New();
    smoother->SetMeshCollection(atlasMeshCollectionPtr);
    smoother->SetSigma(sigma);
    // note: it's very unclear that this actually updates the source mesh collection
    smoother->GetSmoothedMeshCollection();
}

py::array_t<uint16_t> KvlMesh::RasterizeMesh(std::vector<size_t> size, int classNumber) {
    // Some typedefs
    typedef kvl::AtlasMeshAlphaDrawer::ImageType  AlphaImageType;
    typedef AlphaImageType::SizeType  SizeType;
    typedef kvl::AtlasMeshMultiAlphaDrawer::ImageType  MultiAlphasImageType;

    // Retrieve input arguments
    kvl::AtlasMesh::ConstPointer constMesh = static_cast< const kvl::AtlasMesh* >( mesh );
    SizeType imageSize;
    for ( int i = 0; i < 3; i++ )
    {
        imageSize[ i ] = size[i];
    }

//    if ( nrhs > 2 )
//    {
//        double* tmp = mxGetPr( prhs[ 2 ] );
//        classNumber = static_cast< int >( *tmp );
//    }

//std::cout << "mesh: " << mesh.GetPointer() << std::endl;
    //std::cout << "imageSize: " << imageSize << std::endl;
    //std::cout << "classNumber: " << classNumber << std::endl;


    if ( classNumber >= 0 )
    {
//        // Rasterize the specified prior. If the class number is 0, then pre-fill everything
//        // so that parts not overlayed by the mesh are still considered to the background
        kvl::AtlasMeshAlphaDrawer::Pointer  alphaDrawer = kvl::AtlasMeshAlphaDrawer::New();
        alphaDrawer->SetRegions( imageSize );
        alphaDrawer->SetClassNumber( classNumber );
        if ( classNumber == 0 )
        {
            ( const_cast< AlphaImageType* >( alphaDrawer->GetImage() ) )->FillBuffer( 1.0 );
        }
        alphaDrawer->Rasterize( mesh );

//
//
//        // Finally, copy the buffer contents into a Matlab matrix
//        mwSize  dims[ 3 ];
//        for ( int i = 0; i < 3; i++ )
//        {
//            dims[ i ] = imageSize[ i ];
//        }
//        //plhs[ 0 ] = mxCreateNumericArray( 3, dims, mxSINGLE_CLASS, mxREAL );
        std::vector<size_t> shape = {imageSize[ 0 ], imageSize[ 1 ], imageSize[ 2 ]};
//        //float*  data = static_cast< float* >( mxGetData( plhs[ 0 ] ) );
//        plhs[ 0 ] = mxCreateNumericArray( 3, dims, mxUINT16_CLASS, mxREAL );
//        unsigned short*  data = static_cast< unsigned short* >( mxGetData( plhs[ 0 ] ) );
        auto const buffer = new uint16_t[shape[0]*shape[1]*shape[2]];
        auto data = buffer;
//
        itk::ImageRegionConstIterator< AlphaImageType >
                it( alphaDrawer->GetImage(),
                    alphaDrawer->GetImage()->GetBufferedRegion() );
        for ( ;!it.IsAtEnd(); ++it, ++data )
        {
            float  probability = it.Value();
            if ( probability < 0 )
            {
                probability = 0.0f;
            }
            else if ( probability > 1 )
            {
                probability = 1.0f;
            }
            *data = static_cast< unsigned short >( probability * 65535 + .5 );
        }
        return createNumpyArrayFStyle(shape, buffer);

    }
    else
    {
        // Rasterize all priors simultaneously
        const unsigned int  numberOfClasses = mesh->GetPointData()->Begin().Value().m_Alphas.Size();
        //std::cout << "numberOfClasses: " << numberOfClasses << std::endl;

        //std::cout << "Rasterizing mesh..." << std::flush;
        kvl::AtlasMeshMultiAlphaDrawer::Pointer  drawer = kvl::AtlasMeshMultiAlphaDrawer::New();
        drawer->SetRegions( imageSize );
        //std::cout << "here: " << numberOfClasses << std::endl;
        drawer->Rasterize( mesh );
        MultiAlphasImageType::ConstPointer  alphasImage = drawer->GetImage();
        //std::cout << "done" << std::endl;


        // Convert into 4-D Matlab matrix
        //std::cout << "Converting into Matlab format..." << std::flush;

        std::vector<size_t> shape = {imageSize[ 0 ], imageSize[ 1 ], imageSize[ 2 ], numberOfClasses};

        auto const buffer = new uint16_t[shape[0]*shape[1]*shape[2]*shape[3]];
        auto data = buffer;

        for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
        {
            // Loop over all voxels
            itk::ImageRegionConstIterator< MultiAlphasImageType >  alphasIt( alphasImage,
                                                                             alphasImage->GetBufferedRegion() );
            for ( ;!alphasIt.IsAtEnd(); ++alphasIt, ++data )
            {
                float  probability = 0.0f;
                if ( alphasIt.Value().sum() == 0 )
                {
                    // Outside of mesh so not computed. Let's put that to background
                    if ( classNumber == 0 )
                    {
                        //probability = 1.0f;
                    }
                }
                else
                {
                    //*data = alphasIt.Value()[ classNumber ];
                    probability = alphasIt.Value()[ classNumber ];
                    if ( probability < 0 )
                    {
                        probability = 0.0f;
                    }
                    else if ( probability > 1 )
                    {
                        probability = 1.0f;
                    }
                } // End test if outside mesh area

                *data = static_cast< uint16_t >( probability * 65535 + .5 );
            }

            //std::cout << "Done" << std::endl;
        } // End loop over all labels
        //std::cout << "done" << std::endl;
        return createNumpyArrayFStyle(shape, buffer);
    }  // End test if rasterizing one prior or all priors simultaneously

}


py::array KvlMesh::RasterizeValues(std::vector<size_t> size, py::array_t<double, py::array::c_style | py::array::forcecast> values) {

    int nframes = (values.ndim() > 1) ? values.shape(1) : 1;
    
    // get mesh
    kvl::AtlasMesh::ConstPointer constMesh = static_cast< const kvl::AtlasMesh* >( mesh );

    // make value drawer
    typedef kvl::AtlasMeshValueDrawer ValueDrawer;
    ValueDrawer::Pointer valueDrawer = kvl::AtlasMeshValueDrawer::New();

    // set image size
    ValueDrawer::ImageType::SizeType imageSize;
    for (int i = 0; i < 3; i++) imageSize[i] = size[i];
    valueDrawer->SetRegions(imageSize, nframes);
    
    // set values and rasterize
    valueDrawer->SetValues(values.data());
    valueDrawer->Rasterize(mesh);

    // copy to numpy array
    double * const buffer = new double[size[0] * size[1] * size[2] * nframes];
    double * data = buffer;
    for (int frame = 0 ; frame < nframes ; frame++) {
        itk::ImageRegionConstIterator<ValueDrawer::ImageType> it(valueDrawer->GetImage(), valueDrawer->GetImage()->GetBufferedRegion());
        for ( ; !it.IsAtEnd(); ++it, ++data ) *data = it.Value()[frame];
    }

    if (nframes > 1) size.push_back(nframes);
    return createNumpyArrayFStyle(size, buffer);
}


py::array_t<double> KvlMesh::FitAlphas( const py::array_t< uint16_t, py::array::f_style | py::array::forcecast >& probabilityImageBuffer, int EMIterations ) const
{
  
  // Retrieve size of image and number of number of classes
  if ( probabilityImageBuffer.ndim() != 4 ) 
    {
    itkGenericExceptionMacro( "probability image buffer must have four dimensions" );
    }
      
  typedef kvl::AtlasMeshProbabilityImageStatisticsCollector::ProbabilityImageType  ProbabilityImageType; 
  typedef ProbabilityImageType::SizeType  SizeType;
  SizeType  imageSize;
  for ( int i = 0; i < 3; i++ )
    {
    imageSize[ i ] = probabilityImageBuffer.shape( i );
    //std::cout << "imageSize[ i ]: " << imageSize[ i ] << std::endl;
    }
  const int  numberOfClasses = probabilityImageBuffer.shape( 3 );
  std::cout << "imageSize: " << imageSize << std::endl;
  std::cout << "numberOfClasses: " << numberOfClasses << std::endl;


  // Allocate an image of that size
  ProbabilityImageType::Pointer  probabilityImage = ProbabilityImageType::New();
  probabilityImage->SetRegions( imageSize );
  probabilityImage->Allocate();
  ProbabilityImageType::PixelType  emptyEntry( numberOfClasses );
  emptyEntry.Fill( 0.0f );
  probabilityImage->FillBuffer( emptyEntry );
  
  
  // Fill in -- relying on the fact that we've guaranteed a F-style Numpy array input
  auto  data = probabilityImageBuffer.data(); 
  for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
    {
    // Loop over all voxels
    itk::ImageRegionIterator< ProbabilityImageType >  it( probabilityImage,
                                                          probabilityImage->GetBufferedRegion() );
    for ( ;!it.IsAtEnd(); ++it, ++data )
      {
      it.Value()[ classNumber ] = static_cast< float >( *data ) / 65535.0;
      }
    
    }
  std::cout << "Created and filled probabilityImage" << std::endl;  
  

  // EM estimation of probability image in mesh representation. First we create a private mesh
  // to play with - we initialize with a flat prior (alphas)
  
  // Create a flat alpha entry  
  kvl::AtlasAlphasType   flatAlphasEntry( numberOfClasses );
  flatAlphasEntry.Fill( 1.0f / static_cast< float >( numberOfClasses ) );

  // Create a border alphas entry (only first class is possible)
  kvl::AtlasAlphasType   borderAlphasEntry( numberOfClasses );
  borderAlphasEntry.Fill( 0.0f );
  borderAlphasEntry[ 0 ] = 1.0f;
  
  // Initialize point params with flat alphas (unless at border)
  kvl::AtlasMesh::PointDataContainer::Pointer  privateParameters 
                                                 = kvl::AtlasMesh::PointDataContainer::New();
  for ( kvl::AtlasMesh::PointDataContainer::ConstIterator  it = mesh->GetPointData()->Begin();
        it != mesh->GetPointData()->End(); ++it )
    {
    kvl::PointParameters  params = it.Value(); // Copy
    if ( params.m_CanChangeAlphas )
      {
      params.m_Alphas = flatAlphasEntry;
      }
    else
      {
      params.m_Alphas = borderAlphasEntry;
      }
      
    privateParameters->InsertElement( it.Index(), params );     
    }  

  // Create actual private mesh
  kvl::AtlasMesh::Pointer  privateMesh = kvl::AtlasMesh::New();
  kvl::AtlasMesh::Pointer  nonConstMesh = const_cast< kvl::AtlasMesh* >( mesh.GetPointer() );
  privateMesh->SetPoints( nonConstMesh->GetPoints() );
  privateMesh->SetCells( nonConstMesh->GetCells() );
  privateMesh->SetPointData( privateParameters );
  privateMesh->SetCellData( nonConstMesh->GetCellData() );


  // Do the actual EM algorithm using (an updating the alphas in) our private mesh
  for ( int iterationNumber = 0; iterationNumber < EMIterations; iterationNumber++ )
    {
    // E-step: assign voxels to mesh nodes
    kvl::AtlasMeshProbabilityImageStatisticsCollector::Pointer  statisticsCollector = 
                                            kvl::AtlasMeshProbabilityImageStatisticsCollector::New();
    statisticsCollector->SetProbabilityImage( probabilityImage );
    statisticsCollector->Rasterize( privateMesh );
    double  cost = statisticsCollector->GetMinLogLikelihood();
    std::cout << "   EM iteration " << iterationNumber << " -> " << cost << std::endl;
    
    // M-step: normalize class counts in mesh nodes
    kvl::AtlasMesh::PointDataContainer::Iterator  pointParamIt = privateParameters->Begin();
    kvl::AtlasMeshProbabilityImageStatisticsCollector::StatisticsContainerType::ConstIterator  statIt = statisticsCollector->GetLabelStatistics()->Begin();
    for ( ; pointParamIt != privateParameters->End(); ++pointParamIt, ++statIt )
      {
      if ( pointParamIt.Value().m_CanChangeAlphas )
        {
        pointParamIt.Value().m_Alphas = statIt.Value();
        pointParamIt.Value().m_Alphas /= ( statIt.Value().sum() + 1e-12 );
        }
        
      }
  
    } // End loop over EM iteration numbers
   
  // 
  return AlphasToNumpy( privateParameters.GetPointer() ); // Makes a copy
  
}


py::array_t<double> KvlMesh::DrawJacobianDeterminant(std::vector<size_t> size)
{
  // Some typedefs
  typedef kvl::AtlasMeshJacobianDeterminantDrawer::ImageType  JacobianDeterminantImageType;
  typedef JacobianDeterminantImageType::SizeType  SizeType;

  SizeType  imageSize;
  for ( int i = 0; i < 3; i++ )
  {
    imageSize[ i ] = size[ i ];
    //std::cout << "imageSize[ i ]: " << imageSize[ i ] << std::endl;
  }
  std::cout << "imageSize: " << imageSize << std::endl;

  // get mesh
  kvl::AtlasMesh::ConstPointer constMesh = static_cast< const kvl::AtlasMesh* >( mesh );

  // Rasterize
  kvl::AtlasMeshJacobianDeterminantDrawer::Pointer  drawer = kvl::AtlasMeshJacobianDeterminantDrawer::New();
  drawer->SetRegions( imageSize );
  drawer->Rasterize( mesh );

  // Finally, copy the buffer contents into a Numpy array
  auto *outData = new double[imageSize[0] * imageSize[1] * imageSize[2]];
  auto dataIterator = outData;

  itk::ImageRegionConstIterator< JacobianDeterminantImageType > it( drawer->GetImage(), drawer->GetImage()->GetBufferedRegion() );
  for ( ;!it.IsAtEnd(); ++it )
  {
    *dataIterator++ = it.Value();
  }

  return createNumpyArrayFStyle({imageSize[0], imageSize[1], imageSize[2]}, outData);

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

void CreatePointSetFromNumpy(PointSetPointer targetPoints, const py::array_t<double> &source) {

    if (source.ndim() != 2) {
        itkGenericExceptionMacro("point source shape must have two dimensions");
    }
    auto sourceShape = source.shape();
    const unsigned int sourceNumberOfNodes = *sourceShape++;
    const unsigned int sourceXYZDimensions = *sourceShape++;
    if (sourceXYZDimensions != 3) {
        itkGenericExceptionMacro("source points must have 3 coordinates not " << sourceXYZDimensions);
    }
    for (unsigned int pointIndex = 0; pointIndex < sourceNumberOfNodes; ++pointIndex) {
        kvl::AtlasMesh::PointType  point;
        for (int xyzAxisSelector = 0; xyzAxisSelector < XYZ_DIMENSIONS; xyzAxisSelector++) {
            point[xyzAxisSelector] = *source.data(pointIndex, xyzAxisSelector);
        }
        targetPoints->InsertElement( targetPoints->Size(), point );
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

py::array_t<bool> CanMovesToNumpy(PointDataConstPointer pointData ) {
    const unsigned int numberOfNodes = pointData->Size();
    auto *data = new bool[numberOfNodes * 3];
    auto dataIterator = data;
    for (auto srcIterator = pointData->Begin(); srcIterator != pointData->End(); ++srcIterator) {
        *dataIterator++ = srcIterator.Value().m_CanMoveX;
        *dataIterator++ = srcIterator.Value().m_CanMoveY;
        *dataIterator++ = srcIterator.Value().m_CanMoveZ;
    }
    return createNumpyArrayCStyle({numberOfNodes, 3}, data);
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


void CopyNumpyToCanMoves(PointDataPointer destinationPointData, const py::array_t<bool> &source)
{
    if (source.ndim() != 2) {
        itkGenericExceptionMacro("canMove source shape must have two dimensions");
    }
    const unsigned int currentNumberOfNodes = destinationPointData->Size();
    auto sourceShape = source.shape();
    const unsigned int sourceNumberOfNodes = *sourceShape++;
    const unsigned int numberOfDimensions = *sourceShape++;
    if (sourceNumberOfNodes != currentNumberOfNodes) {
        itkGenericExceptionMacro(
                "source data point count of "
                        << sourceNumberOfNodes
                        << " not equal to mesh point count of "
                        << currentNumberOfNodes
        );
    }
    if (numberOfDimensions !=3 ) {
        itkGenericExceptionMacro("source data have three directions not " << numberOfDimensions);
    }
    unsigned int pointIndex = 0;
    for (auto destIterator = destinationPointData->Begin(); destIterator != destinationPointData->End(); ++destIterator, ++pointIndex) {
        destIterator.Value().m_CanMoveX = *source.data(pointIndex, 0);
        destIterator.Value().m_CanMoveY = *source.data(pointIndex, 1);
        destIterator.Value().m_CanMoveZ = *source.data(pointIndex, 2);
    }
}


KvlMesh* KvlMesh::GetSubmesh( py::array_t<bool>& mask ) {
  
  if (mask.ndim() < 1 ) 
    {
    itkGenericExceptionMacro("mask shape must have at least one dimension");
    }
  const unsigned int  numberOfMeshNodes = mesh->GetPoints()->Size();
  const unsigned int  maskNumberOfNodes = mask.shape( 0 );
  if ( numberOfMeshNodes != maskNumberOfNodes )
    {
    itkGenericExceptionMacro( "Numbe of nodes in mask not compatible with number of nodes in mesh (" 
                              << maskNumberOfNodes << " vs. " << numberOfMeshNodes << ")" );
    }
    
  // Collect a list of points that the mask indicates should be included    
  std::set< kvl::AtlasMesh::PointIdentifier >  initialPointIds;
  int  pointNumber = 0;
  for ( auto it = mesh->GetPoints()->Begin(); it != mesh->GetPoints()->End(); 
        ++it, ++pointNumber ) 
    {
    if ( *( mask.data( pointNumber ) ) )
      {
      //const int  counter = pointIdLookupTable.size();
      //pointIdLookupTable[ it.Index() ] = counter;
      initialPointIds.insert( it.Index() );
      }
    }
    
  
  // Build a list of tetrahedra that are attached to these nodes. At the same
  // time, create a simple lookup table to converts the original cellIds into
  // new cellIds (without gaps) in our new mesh 
  std::map< kvl::AtlasMesh::CellIdentifier, int >  tetIdLookupTable;
  mesh->BuildCellLinks();
  const kvl::AtlasMesh::CellLinksContainer::ConstPointer  cellLinks = mesh->GetCellLinks();
  for ( auto pointIt = initialPointIds.begin(); pointIt != initialPointIds.end(); ++pointIt )
    {
    // Loop over all the tetrahedra attached to this point  
    const std::set< kvl::AtlasMesh::CellIdentifier >&  cells 
                                                       = cellLinks->ElementAt( *pointIt );
    for ( auto cellIt = cells.begin(); cellIt != cells.end(); ++cellIt )
      {
      //  
      const kvl::AtlasMesh::CellIdentifier  cellId = *cellIt;
      const kvl::AtlasMesh::CellType*  cell = mesh->GetCells()->ElementAt( cellId );
      if ( cell->GetType() != kvl::AtlasMesh::CellType::TETRAHEDRON_CELL )
        {
        //  
        continue;  
        }
        
      // OK, we're having a tetrahedron. Add it to our list if not already there
      if ( tetIdLookupTable.find( cellId ) == tetIdLookupTable.end() )
        {
        const int  counter = tetIdLookupTable.size();
        tetIdLookupTable[ cellId ] = counter;
        }
      
      } // End loop over cells/tets
  
    } // End loop over points
  

  // The tetrahedra will have extra points (not in the original mask) attached to them.
  // These should also be part of our new mesh.
  std::set< kvl::AtlasMesh::PointIdentifier >  extraPointIds;
  for ( auto tetIt = tetIdLookupTable.begin(); tetIt != tetIdLookupTable.end(); ++tetIt )
    {
    //  
    //const kvl::AtlasMesh::CellType&  tet = mesh->GetCells()->ElementAt( tetIt->first );
    const kvl::AtlasMesh::CellType*  tet = mesh->GetCells()->ElementAt( tetIt->first );
  
  
    // Loop over all points in tet
    for ( auto pit = tet->PointIdsBegin(); pit != tet->PointIdsEnd(); ++pit )
      {
      const kvl::AtlasMesh::PointIdentifier  pointId = *pit;
        
      if ( initialPointIds.find( pointId ) == initialPointIds.end() )
        {
        // Found a new, extra point  
        extraPointIds.insert( pointId );
        }
      } // End loop over all points in tet
      
    } // End loop over all tets  
    
    
  // Build a look-up table converting the pointIds of all the relevant points in the 
  // original mesh into a new pointId numbering system (without gaps) in our new mesh.
  std::map< kvl::AtlasMesh::PointIdentifier, int >  pointIdLookupTable;
  std::set< kvl::AtlasMesh::PointIdentifier >  mergedPointIds( initialPointIds );
  mergedPointIds.insert( extraPointIds.begin(), extraPointIds.end() );
  for ( auto it = mergedPointIds.begin(); it != mergedPointIds.end(); ++it )
    {
    const int  counter = pointIdLookupTable.size();
    pointIdLookupTable[ *it ] = counter;
    }
  
  
  // Copy points and pointData   
  kvl::AtlasMesh::PointsContainer::Pointer  subPoints =  kvl::AtlasMesh::PointsContainer::New();
  kvl::AtlasMesh::PointDataContainer::Pointer  subPointData 
                                                      = kvl::AtlasMesh::PointDataContainer::New();
  for ( auto pointIt = pointIdLookupTable.begin(); pointIt != pointIdLookupTable.end(); ++pointIt )
    {
    subPoints->InsertElement( pointIt->second, mesh->GetPoints()->ElementAt( pointIt->first ) );
    subPointData->InsertElement( pointIt->second, mesh->GetPointData()->ElementAt( pointIt->first ) );
    }
    
    
  // Same for cells and cellData
  kvl::AtlasMesh::CellsContainer::Pointer  subCells = kvl::AtlasMesh::CellsContainer::New();
  kvl::AtlasMesh::CellDataContainer::Pointer  subCellData = kvl::AtlasMesh::CellDataContainer::New();
  for ( auto tetIt = tetIdLookupTable.begin(); tetIt != tetIdLookupTable.end(); ++tetIt )
    {
    // Create new tet cell  
    typedef itk::TetrahedronCell< kvl::AtlasMesh::CellType >  TetrahedronCell;
    kvl::AtlasMesh::CellAutoPointer  newCell;
    newCell.TakeOwnership( new TetrahedronCell );
    const kvl::AtlasMesh::CellType*  tet = mesh->GetCells()->ElementAt( tetIt->first );

    int  localId = 0;
    for ( auto pit = tet->PointIdsBegin(); pit != tet->PointIdsEnd(); ++pit, ++localId )
      {
      const kvl::AtlasMesh::PointIdentifier  pointId = *pit;
      newCell->SetPointId( localId, pointIdLookupTable[ pointId ] );
      }

    // Insert new cell  
    subCells->InsertElement( tetIt->second, newCell.ReleaseOwnership() );
  
    // Also copy cell data
    subCellData->InsertElement( tetIt->second, mesh->GetCellData()->ElementAt( tetIt->first ) );

    } // End loop over tets
  
  
  // Last but not least: the extra points introduced because they are part of boundary tets
  // need to be set to immobile. Finding them is easy because their pointNumber exceeds the
  // original number of points
  for ( auto pointIt = extraPointIds.begin(); pointIt != extraPointIds.end(); ++pointIt )
    {
    const int  pointNumber = pointIdLookupTable[ *pointIt ];  
    subPointData->ElementAt( pointNumber ).m_CanMoveX = false;
    subPointData->ElementAt( pointNumber ).m_CanMoveY = false;
    subPointData->ElementAt( pointNumber ).m_CanMoveZ = false;
    }  
  
  
  // The original mask will not yet have the extra points in it -- add back as
  // feedback to the user
  auto  x = mask.mutable_unchecked();
  pointNumber = 0;
  for ( auto it = mesh->GetPoints()->Begin(); it != mesh->GetPoints()->End(); 
        ++it, ++pointNumber ) 
    {
    if ( pointIdLookupTable.find( it.Index() ) != pointIdLookupTable.end() )
      {
      x( pointNumber ) = true;
      }
    }
  
  
  // Finally, create a mesh
  kvl::AtlasMesh::Pointer  subMesh = kvl::AtlasMesh::New();
  subMesh->SetPoints( subPoints );
  subMesh->SetPointData( subPointData );
  subMesh->SetCells( subCells );
  subMesh->SetCellData( subCellData );


  //
  kvl::AtlasMesh::ConstPointer  tmp = subMesh.GetPointer();
  return new KvlMesh( tmp );
  
}

