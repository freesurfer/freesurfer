#include "kvlAtlasMesh.h"
#include "kvlAtlasMeshCollection.h"
#include "kvlAtlasMeshDeformationLBFGSOptimizer.h"
#include "kvlAtlasMeshToIntensityImageCostAndGradientCalculator.h"
#include "itkCommand.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkMGHImageIOFactory.h"
#include "kvlCroppedImageReader.h"

// Usage: 
//   testdeformMeshPython <maxIter> <numThreads> <input-mesh-collection>  <output-deformed-mesh-collection> <input-mgz>
// examples:
//   testdeformMeshPython 30 1 atlas_level1_deformmesh.txt.gz  deformedmesh_out.txt.gz inp_image_deformmesh.mgz
int main( int argc, char** argv )
{
  // Add support for MGH file format to ITK. An alternative way to add this by default would be
  // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
  itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );

  // Retrieve the input parameters
  std::ostringstream  inputParserStream;
  for ( int argumentNumber = 1; argumentNumber < 5; argumentNumber++ ) 
    {
    inputParserStream << argv[ argumentNumber ] << " ";
    }
  std::istringstream  inputStream( inputParserStream.str().c_str() );
  int maxIteration;
  int numThreads;
  std::string  meshCollectionFile;
  std::string  deformedMeshCollectionFile;
  inputStream >> maxIteration >> numThreads >> meshCollectionFile >> deformedMeshCollectionFile;

  std::cout << "deformMesh Command line params:" << std::endl;
  std::cout << "  maxIteration:                         " << maxIteration << std::endl;
  std::cout << "  numThreads:                           " << numThreads   << std::endl;
  std::cout << "  meshCollectionFile:                   " << meshCollectionFile << std::endl;
  std::cout << "  output deformedMeshCollectionFile:    " << deformedMeshCollectionFile << std::endl;

  // set ITK number of threads  
  std::cout << "[DEBUG] itk::MultiThreader::SetGlobalDefaultNumberOfThreads(" << numThreads << ")" << std::endl;
  itk::MultiThreader::SetGlobalDefaultNumberOfThreads(numThreads);

  // read mesh collection
  kvl::AtlasMeshCollection::Pointer  meshStartCollection = 0;
  if ( itksys::SystemTools::FileExists( meshCollectionFile.c_str(), true ) )
    {
    meshStartCollection = kvl::AtlasMeshCollection::New();
    if ( !meshStartCollection->Read( meshCollectionFile.c_str() ) )
      {
      std::cerr << "Couldn't read mesh from file " << meshCollectionFile << std::endl;
      return -1;
      } 
    else
      {
      std::cout << "[DEBUG] " << meshCollectionFile << " Found!" << std::endl;
      std::cout << "meshStartCollection found; reading from: " << meshCollectionFile << std::endl;
      }
    }

  typedef kvl::CroppedImageReader::TransformType TransformType;
  typedef TransformType::Pointer TransformPointer;

  /**************************/
  /* read images            */
  /**************************/
  const std::string  imageFileName = "/autofs/space/curv_001/users/yujing/input.mgz";
  //typedef kvl::CompressionLookupTable::ImageType  LabelImageType;     // typedef itk::Image< unsigned short, 3 >  ImageType;
  typedef itk::Image< float, 3 > ImageType;
  std::vector< ImageType::ConstPointer >  images;
  for ( int argumentNumber = 5; argumentNumber < argc; argumentNumber++ )
  {
    std::string imageFileName = argv[ argumentNumber ];

#if 0
    // Read the image  (KvlImage::KvlImage(const std::string &imageFileName))
    kvl::CroppedImageReader::Pointer reader = kvl::CroppedImageReader::New();
    reader->Read( imageFileName.c_str() );

    // Convert the image to float
    typedef itk::CastImageFilter< kvl::CroppedImageReader::ImageType, ImageType > CasterType;
    CasterType::Pointer  caster = CasterType::New();
    caster->SetInput( reader->GetImage() );
    caster->Update();

    // Store the image and transform in persistent memory
    ImageType::ConstPointer  image = caster->GetOutput();
    TransformPointer transform = TransformType::New();
    reader->GetWorldToImageTransform()->GetInverse( transform );
    std::cout << "Read image: " << imageFileName << std::endl;

#else
    std::cout << "Reading input image: " << argv[ argumentNumber ] << std::endl;

    // Read the input image
    typedef itk::ImageFileReader< ImageType >  ReaderType;
    ReaderType::Pointer  reader = ReaderType::New();
    reader->SetFileName( imageFileName );
    reader->Update();
    ImageType::ConstPointer  image = reader->GetOutput();

#if 1
    // Over-ride the spacing and origin since at this point we can't deal with that
    const double spacing[] = { 1, 1, 1 };
    const double origin[] = { 0, 0, 0 };
    const_cast< ImageType* >( image.GetPointer() )->SetSpacing( spacing );
    const_cast< ImageType* >( image.GetPointer() )->SetOrigin( origin );
#endif
#endif

    // Remember this image
    images.push_back( image );
  }

  /*******************************************************************************/
  /* set up optimizer                                                            */
  /*                                                                             */
  /*      optimizerType = 'L-BFGS'                                               */
  /*      optimizationParameters = {                                             */
  /*          'Verbose': False,                                                  */
  /*          'MaximalDeformationStopCriterion': 0.001,  # measured in pixels,   */
  /*          'LineSearchMaximalDeformationIntervalStopCriterion': 0.001,        */
  /*          'MaximumNumberOfIterations': 20,                                   */
  /*          'BFGS-MaximumMemoryLength': 12                                     */
  /*******************************************************************************/
  bool  verbose = false;
  float maximalDeformationStopCriterion = 0.001;
  float lineSearchMaximalDeformationIntervalStopCriterion = 0.001;
  int   maximumNumberOfIterations = maxIteration;  //20;
  int   maximumMemoryLength = 12;
  kvl::AtlasMeshDeformationLBFGSOptimizer::Pointer  optimizer
    = kvl::AtlasMeshDeformationLBFGSOptimizer::New();
  optimizer->SetVerbose( true );
  optimizer->SetMaximalDeformationStopCriterion( maximalDeformationStopCriterion );
  optimizer->SetMaximumNumberOfIterations( maximumNumberOfIterations );
  optimizer->SetLineSearchMaximalDeformationIntervalStopCriterion( lineSearchMaximalDeformationIntervalStopCriterion );
  optimizer->SetMaximumMemoryLength( maximumMemoryLength );

  /***********************************************************************/
  /*  Set up gradient calculator                                         */
  /*                                                                     */
  /*          typeName='AtlasMeshToIntensityImage',                      */
  /*          images=images,                                             */
  /*          boundaryCondition='Sliding',                               */
  /*          transform=transform,                                       */
  /*          means=means,                                               */
  /*          variances=variances,                                       */
  /*          mixtureWeights=mixtureWeights,                             */
  /*          numberOfGaussiansPerClass=numberOfGaussiansPerClass)       */
  /***********************************************************************/
  kvl::AtlasMeshToIntensityImageCostAndGradientCalculator::Pointer calculator
    = kvl::AtlasMeshToIntensityImageCostAndGradientCalculator::New();

  // hard-coded calculator parameters
  const int  numberOfGaussians = 9;
  const int  numberOfContrasts = 1;
  std::vector< vnl_vector< double > >  means_converted;
  vnl_vector< double > v( numberOfContrasts, 0.0f );
  v[0] = 2.03099;  means_converted.push_back(v);
  v[0] = 5.01833;  means_converted.push_back(v);
  v[0] = 4.63988;  means_converted.push_back(v);
  v[0] = 4.04694;  means_converted.push_back(v);
  v[0] = 3.82442;  means_converted.push_back(v);
  v[0] = 4.90528;  means_converted.push_back(v);
  v[0] = 4.08647;  means_converted.push_back(v);
  v[0] = 4.63836;  means_converted.push_back(v);
  v[0] = 4.52726;  means_converted.push_back(v);
#if 0
  std::cout << "[DEBUG] means_converted.size=" << means_converted.size() << std::endl;
  std::cout << "[DEBUG] means_converted[ 0 ].size()=" << means_converted[ 0 ].size() << std::endl;
  for ( int gaussianNumber = 0; gaussianNumber < numberOfGaussians; gaussianNumber++ )
  {
    std::cout << "[DEBUG] means_converted[" << gaussianNumber << "]=" << means_converted.at(gaussianNumber) << std::endl;
  }
#endif

  std::vector< vnl_matrix< double > >  variances_converted;
  vnl_matrix< double >  variance( numberOfContrasts, numberOfContrasts, 0.0f );
  variance[0][0] = 0.645518;    variances_converted.push_back( variance );
  variance[0][0] = 0.0133676;   variances_converted.push_back( variance );
  variance[0][0] = 0.0350092;   variances_converted.push_back( variance );
  variance[0][0] = 0.188636;    variances_converted.push_back( variance );
  variance[0][0] = 1.01978;     variances_converted.push_back( variance );
  variance[0][0] = 0.00798006;  variances_converted.push_back( variance );
  variance[0][0] = 0.4097;      variances_converted.push_back( variance );
  variance[0][0] = 0.164445;    variances_converted.push_back( variance );
  variance[0][0] = 0.0607125;   variances_converted.push_back( variance );
#if 0
  std::cout << "[DEBUG] variances_converted[ 0 ].rows()=" << variances_converted[ 0 ].rows() << std::endl;
  for ( unsigned int gaussianNumber = 0; gaussianNumber < numberOfGaussians; gaussianNumber++ ) 
  {
    std::cout << "[DEBUG] variances_converted[" << gaussianNumber << "]=" << variances_converted.at(gaussianNumber) << std::endl;
  }
#endif

  std::vector< double >  mixtureWeights_converted;
  double mixtureWeight;
  mixtureWeight = 1;  mixtureWeights_converted.push_back(mixtureWeight);
  mixtureWeight = 1;  mixtureWeights_converted.push_back(mixtureWeight);
  mixtureWeight = 1;  mixtureWeights_converted.push_back(mixtureWeight);
  mixtureWeight = 1;  mixtureWeights_converted.push_back(mixtureWeight);
  mixtureWeight = 1;  mixtureWeights_converted.push_back(mixtureWeight);
  mixtureWeight = 1;  mixtureWeights_converted.push_back(mixtureWeight);
  mixtureWeight = 0.186982;  mixtureWeights_converted.push_back(mixtureWeight);
  mixtureWeight = 0.367688;  mixtureWeights_converted.push_back(mixtureWeight);
  mixtureWeight = 0.44533;  mixtureWeights_converted.push_back(mixtureWeight);
#if 0
  std::cout << "[DEBUG] mixtureWeights_converted.size()=" << mixtureWeights_converted.size() << ",numberOfGaussians=" << numberOfGaussians << std::endl;
  for ( int gaussianNumber = 0; gaussianNumber < numberOfGaussians; gaussianNumber++ )
  {
    std::cout << "[DEBUG] mixtureWeights_converted[" << gaussianNumber << "] = " << mixtureWeights_converted.at(gaussianNumber) << std::endl;
  }
#endif

  const int  numberOfClasses = 7;
  std::vector< int >  numberOfGaussiansPerClass_converted;
  numberOfGaussiansPerClass_converted.push_back(1);
  numberOfGaussiansPerClass_converted.push_back(1);
  numberOfGaussiansPerClass_converted.push_back(1);
  numberOfGaussiansPerClass_converted.push_back(1);
  numberOfGaussiansPerClass_converted.push_back(1);
  numberOfGaussiansPerClass_converted.push_back(1);
  numberOfGaussiansPerClass_converted.push_back(3);
#if 0
  int  sum = 0;
  for ( int i = 0; i < numberOfGaussiansPerClass_converted.size(); i++ )
  {
    sum += numberOfGaussiansPerClass_converted[ i ];
  }
  std::cout << "[DEBUG] numberOfGaussiansPerClass_converted.size()=" << numberOfGaussiansPerClass_converted.size() << ",sum=" << sum << ",numberOfClasses=" << numberOfClasses << std::endl;
 
  for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
  {
    std::cout << "[DEBUG] numberOfGaussiansPerClass_converted[" << classNumber << "] = " << numberOfGaussiansPerClass_converted.at(classNumber) << std::endl;
  }
#endif

  calculator->SetImages( images );
  calculator->SetParameters( means_converted, variances_converted, mixtureWeights_converted, numberOfGaussiansPerClass_converted );
  calculator->SetBoundaryCondition( kvl::AtlasMeshPositionCostAndGradientCalculator::SLIDING );

  // set optimizer mesh
  kvl::AtlasMesh::ConstPointer constMesh = meshStartCollection->GetReferenceMesh();
  kvl::AtlasMesh::Pointer mutableMesh = const_cast< kvl::AtlasMesh* >( constMesh.GetPointer() );
  // Pass the mesh and calculator to it
  optimizer->SetMesh( mutableMesh );
  optimizer->SetCostAndGradientCalculator( calculator );

  std::cout << "[DEBUG] Done setup AtlasMeshDeformationLBFGSOptimizer and AtlasMeshToIntensityImageCostAndGradientCalculator" << std::endl;
  std::cout << std::endl;
  int iter = 0, totalRasterizeCalls = 0, totalTet = 0, totalZeroVoxelTet = 0, totalVoxel = 0, totalVoxelInTetrahedron = 0;
  while (true)      //while (iter < maxIteration)
  {
#ifdef GEMS_DEBUG_RASTERIZE_VOXEL_COUNT
    calculator->m_Iterations++; 
#endif

    const double  maximalDeformation = optimizer->Step();
    const double  minLogLikelihoodTimesPrior = optimizer->GetMinLogLikelihoodTimesPrior();
    std::cout << "[DEBUG] iter. #" << iter+1 << ", maximalDeformation = " << maximalDeformation << ", minLogLikelihoodTimesPrior = " << minLogLikelihoodTimesPrior << std::endl;
    std::cout << std::endl;

#ifdef GEMS_DEBUG_RASTERIZE_VOXEL_COUNT
    printf("[DEBUG] %3d Rasterize() calls, total tetrahedron = %8d, zero voxel tetrahedron     = %8d, %6.2f\%\n"
           "                               total voxel       = %8d, total voxel in tetrahedron = %8d, %6.2f\%\n\n", 
           calculator->m_RasterizeCalls, calculator->m_tetrahedronCnt, calculator->m_zeroVoxel_tetrahedronCnt, 
           (calculator->m_tetrahedronCnt == 0) ? 0.0 : (double)calculator->m_zeroVoxel_tetrahedronCnt/calculator->m_tetrahedronCnt*100,
           calculator->m_totalVoxel, calculator->m_totalVoxelInTetrahedron, 
           (calculator->m_totalVoxel == 0) ? 0.0 : (double)calculator->m_totalVoxelInTetrahedron/calculator->m_totalVoxel*100);

    totalRasterizeCalls += calculator->m_RasterizeCalls;
    calculator->m_RasterizeCalls = 0; 

    totalTet += calculator->m_tetrahedronCnt; totalZeroVoxelTet += calculator->m_zeroVoxel_tetrahedronCnt;
    calculator->m_tetrahedronCnt = 0; calculator->m_zeroVoxel_tetrahedronCnt = 0;

    totalVoxel += calculator->m_totalVoxel; totalVoxelInTetrahedron += calculator->m_totalVoxelInTetrahedron;
    calculator->m_totalVoxel = 0, calculator->m_totalVoxelInTetrahedron = 0;
#endif

    if (maximalDeformation == 0)
      break;

    iter++;
  }

#ifdef GEMS_DEBUG_RASTERIZE_VOXEL_COUNT
  printf("[DEBUG] Summary: Total %3d Rasterize() calls, total tetrahedron = %12d, zero voxel tetrahedron     = %12d, %6.2f\%\n"
         "                                              total voxel       = %12d, total voxel in tetrahedron = %12d, %6.2f\%\n\n", 
         totalRasterizeCalls, totalTet, totalZeroVoxelTet, 
         (totalTet == 0) ? 0.0 : (double)totalZeroVoxelTet/totalTet*100,
         totalVoxel, totalVoxelInTetrahedron, 
         (totalVoxel == 0) ? 0.0 : (double)totalVoxelInTetrahedron/totalVoxel*100);
#endif

  // output the deformed mesh 
  kvl::AtlasMeshCollection::Pointer  deformedMeshCollection = kvl::AtlasMeshCollection::New();
  deformedMeshCollection->GenerateFromSingleMesh(mutableMesh, 1, 0.1);
  deformedMeshCollection->Write(deformedMeshCollectionFile.c_str());

  return 0;
};

