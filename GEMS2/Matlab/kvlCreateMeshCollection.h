#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "kvlAtlasMeshCollection.h"
#include "kvlAtlasMeshProbabilityImageStatisticsCollector.h"


namespace kvl
{

  
class CreateMeshCollection : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef CreateMeshCollection         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( CreateMeshCollection, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
              
    // meshCollection = kvlCreateMeshCollection( multiAlphaImageBuffer, meshSize )
  
    // Make sure input arguments are correct
    if ( ( nrhs != 2 ) || ( nlhs != 1 ) || 
         ( mxGetNumberOfDimensions( prhs[ 0 ] ) != 4 ) ||
         ( !mxIsUint16( prhs[ 0 ] ) ) )
      {
      mexErrMsgTxt( "Incorrect arguments" );
      }
      
      
      
    // Determine the size of the image to be created  
    typedef AtlasMeshProbabilityImageStatisticsCollector::ProbabilityImageType  ProbabilityImageType; 
    typedef ProbabilityImageType::SizeType  SizeType;
    SizeType  imageSize;
    for ( int i = 0; i < 3; i++ )
      {
      imageSize[ i ] = mxGetDimensions( prhs[ 0 ] )[ i ];
      std::cout << "imageSize[ i ]: " << imageSize[ i ] << std::endl;
      }
 
    // Allocate an image of that size
    const int  numberOfClasses = mxGetDimensions( prhs[ 0 ] )[ 3 ];
    std::cout << "numberOfClasses: " << numberOfClasses << std::endl;
    ProbabilityImageType::Pointer  probabilityImage = ProbabilityImageType::New();
    probabilityImage->SetRegions( imageSize );
    probabilityImage->Allocate();
    ProbabilityImageType::PixelType  emptyEntry( numberOfClasses );
    emptyEntry.Fill( 0.0f );
    probabilityImage->FillBuffer( emptyEntry );


    // Fill in  
    unsigned short*  data = static_cast< unsigned short* >( mxGetData( prhs[ 0 ] ) ); 
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
      
      
    //
    unsigned int  meshSize[ 3 ];
    double*  tmp = mxGetPr( prhs[ 1 ] );
    for ( int i = 0; i < 3; i++, tmp++ )
      {
      meshSize[ i ] = static_cast< int >( *tmp );
      std::cout << "meshSize[ " << i << " ]: " << meshSize[ i ] << std::endl;
      }
    unsigned int  domainSize[ 3 ];
    domainSize[ 0 ] = imageSize[ 0 ];
    domainSize[ 1 ] = imageSize[ 1 ];
    domainSize[ 2 ] = imageSize[ 2 ];

    const float  initialStiffness = 0.1;
    const unsigned int  numberOfMeshes = 1;
    AtlasMeshCollection::Pointer  meshCollection = AtlasMeshCollection::New();
    meshCollection->Construct( meshSize, domainSize, initialStiffness, 
                               numberOfClasses, numberOfMeshes );
    
    std::cout << "Created a mesh collection" << std::endl;
    
    
      
    // EM estimation of the mesh node alphas that best fit the probabilityImage
    //meshCollection->FlattenAlphas(); // Initialization to flat alphas
    std::cout << "Estimating mesh collection alphas using EM" << std::endl;
    for ( int iterationNumber = 0; iterationNumber < 10; iterationNumber++ )
      {
      // E-step
      AtlasMeshProbabilityImageStatisticsCollector::Pointer  statisticsCollector = 
                                              AtlasMeshProbabilityImageStatisticsCollector::New();
      statisticsCollector->SetProbabilityImage( probabilityImage );
      statisticsCollector->Rasterize( meshCollection->GetReferenceMesh() );
      const double  cost = statisticsCollector->GetMinLogLikelihood();
      std::cout << "   EM iteration " << iterationNumber << " -> " << cost << std::endl;
      
      // M-step: Normalize and assign to the meshCollection's alpha vectors
      AtlasMesh::PointDataContainer::Iterator  pointParamIt = meshCollection->GetPointParameters()->Begin();
      AtlasMeshProbabilityImageStatisticsCollector::StatisticsContainerType::ConstIterator  
                                               statIt = statisticsCollector->GetLabelStatistics()->Begin();
      for ( ; pointParamIt != meshCollection->GetPointParameters()->End(); ++pointParamIt, ++statIt )
        {
        if ( pointParamIt.Value().m_CanChangeAlphas )
          {
          pointParamIt.Value().m_Alphas = statIt.Value();
          pointParamIt.Value().m_Alphas /= ( statIt.Value().sum() + 1e-12 );
          }
          
        }
    
      } // End loop over EM iteration numbers
    std::cout << "Done!" << std::endl;    
    
    
    
    // Store the created mesh collection in persistent memory
    const int  meshCollectionHandle = kvl::MatlabObjectArray::GetInstance()->AddObject( meshCollection );
     
    // Return the handle to Matlab
    mwSize  dims[ 1 ];
    dims[ 0 ] = 1;
    plhs[ 0 ] = mxCreateNumericArray( 1, dims, mxINT64_CLASS, mxREAL );
    *( static_cast< int* >( mxGetData( plhs[ 0 ] ) ) ) = meshCollectionHandle;


    }  
  
protected:
  CreateMeshCollection() {};
  virtual ~CreateMeshCollection() {};


  CreateMeshCollection(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl


