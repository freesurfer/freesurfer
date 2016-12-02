#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlAtlasMeshToIntensityImageGradientCalculator.h"


namespace kvl
{

class EvaluateMeshPosition : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef EvaluateMeshPosition         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( EvaluateMeshPosition, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
              
    // [ cost gradient ] = kvlEvaluateMeshPosition( mesh, images, transform, means, precisions ) 
  
    // Make sure input arguments are correct
    if ( ( nrhs < 5 ) || 
         !mxIsInt64( prhs[ 0 ] ) || 
         !mxIsInt64( prhs[ 1 ] ) || 
         !mxIsInt64( prhs[ 2 ] ) ||
         !mxIsDouble( prhs[ 3 ] ) ||
         !mxIsDouble( prhs[ 4 ] ) )
      {
      mexErrMsgTxt( "Incorrect arguments" );
      }
      
    // Retrieve input mesh
    const int meshHandle = *( static_cast< int* >( mxGetData( prhs[ 0 ] ) ) );
    itk::Object::ConstPointer object = kvl::MatlabObjectArray::GetInstance()->GetObject( meshHandle );
    if ( typeid( *object ) != typeid( kvl::AtlasMesh ) )
      {
      mexErrMsgTxt( "mesh doesn't refer to the correct ITK object type" );
      }
    kvl::AtlasMesh::ConstPointer constMesh = static_cast< const kvl::AtlasMesh* >( object.GetPointer() );
    kvl::AtlasMesh::Pointer mesh = const_cast< kvl::AtlasMesh* >( constMesh.GetPointer() );

    // Retrieve input image(s)
    //typedef itk::Image< unsigned short, 3 >  ImageType;
    typedef AtlasMeshDeformationConjugateGradientOptimizer::ImageType  ImageType;

    const int  N = mxGetN( prhs[ 1 ] );
    const int  M = mxGetM( prhs[ 1 ] );
    int numberOfImages = 0;
    
    if(N<M)
    {
      numberOfImages = M;
    }
    else
    {
      numberOfImages = N;
    }

    mexPrintf("numberOfImages = %d\n",numberOfImages);
    std::vector<ImageType::Pointer> images;
    uint64_T *  imagesHandle =  static_cast< uint64_T * >( mxGetData( prhs[ 1 ] ) );

    for(unsigned int nima = 0; nima < numberOfImages; nima++, imagesHandle++)
       {
	 
	 const int handle = *(imagesHandle);
	 std::cout<<"Image: "<<handle<<std::endl;
	 itk::Object::ConstPointer object = kvl::MatlabObjectArray::GetInstance()->GetObject( handle );
	 if ( typeid(*(object)  ) != typeid( ImageType ) )
	   {
	     mexErrMsgTxt( "image doesn't refer to the correct ITK object type" );
	   }
	 ImageType::ConstPointer constImage = static_cast< const ImageType* >( object.GetPointer() );
         ImageType::Pointer image = const_cast< ImageType* >( constImage.GetPointer() );
         images.push_back(image); 

       }
    
    
    // Retrieve transform
    typedef CroppedImageReader::TransformType  TransformType;
    const int transformHandle = *( static_cast< int* >( mxGetData( prhs[ 2 ] ) ) );
    object = kvl::MatlabObjectArray::GetInstance()->GetObject( transformHandle );
    if ( typeid( *object ) != typeid( TransformType ) )
      {
      mexErrMsgTxt( "transform doesn't refer to the correct ITK object type" );
      }
    TransformType::ConstPointer  constTransform = static_cast< const TransformType* >( object.GetPointer() );
    TransformType::Pointer  transform = const_cast< TransformType* >( constTransform.GetPointer() );

   
    // Retrieve means and variances
    const int  numberOfClasses = mxGetDimensions( prhs[ 3 ] )[ 0 ];
    std::vector< vnl_vector< float > > means;
    std::vector< vnl_matrix< float > > precisions;
    vnl_vector< float >  mean ( numberOfImages, 0.0f );
    vnl_matrix< float >  precision( numberOfImages, numberOfImages, 0.0f);

    for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
	for ( int nima = 0; nima < numberOfImages; nima++ )
	  {
	    mean[ nima ] = (mxGetPr( prhs[ 3 ] ))[ classNumber + numberOfClasses*nima ];
	   
	  }
	means.push_back(mean);
      }
    
    //Does not really matter which way you read these in, because the precisions are symmetric matrices
    //transpose wont do any harm. 
    for ( unsigned int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
    	for ( unsigned int row = 0; row < numberOfImages; row++ )
    	  {
    	    for ( unsigned int col = 0; col < numberOfImages; col++ )
    	      {
    		precision[ row ][ col ] = mxGetPr( prhs[ 4 ] )[ row + numberOfImages*(col + numberOfImages*classNumber) ];
    	      }
   	  }
    	precisions.push_back(precision);
      }
  
       

    // Show what we have so far
    //std::cout << "mesh: " << mesh.GetPointer() << std::endl;
    //std::cout << "image: " << image.GetPointer() << std::endl;
    //std::cout << "transform: " << transform.GetPointer() << std::endl;
    //std::cout << "means: " << means << std::endl;
    //std::cout << "variances: " << variances << std::endl;
    
    
    // Set up the gradient calculator
    AtlasMeshToIntensityImageGradientCalculator::Pointer  gradientCalculator 
                              = AtlasMeshToIntensityImageGradientCalculator::New();
    typedef AtlasMeshToIntensityImageGradientCalculator::LabelImageType  DummyTemplateImageType;
    DummyTemplateImageType::Pointer  dummyTemplateImage = DummyTemplateImageType::New();
    dummyTemplateImage->SetRegions( images[0]->GetBufferedRegion() );
    //dummyTemplateImage->Allocate();
         
    gradientCalculator->SetLabelImage( dummyTemplateImage );
    gradientCalculator->SetImages( images );
    gradientCalculator->SetMeans( means );
    gradientCalculator->SetPrecisions( precisions );
    gradientCalculator->SetMeshToImageTransform( transform );
    
    // Let the beast go
    gradientCalculator->Rasterize( mesh );

    // Retrieve the result
    const double  cost = gradientCalculator->GetMinLogLikelihoodTimesPrior();
    AtlasPositionGradientContainerType::ConstPointer  
          gradient = gradientCalculator->GetPositionGradient().GetPointer();
    
    
    // Return the cost and gradient to Matlab
    plhs[ 0 ] = mxCreateDoubleScalar( cost );

    const int  numberOfNodes = gradient->Size();
    //std::cout << "numberOfNodes :" << numberOfNodes << std::endl;
    mwSize  dims[ 2 ];
    dims[ 0 ] = numberOfNodes;
    dims[ 1 ] = 3;
    plhs[ 1 ] = mxCreateNumericArray( 2, dims, mxDOUBLE_CLASS, mxREAL );
    double*  data = static_cast< double* >( mxGetData( plhs[ 1 ] ) ); 

    for ( AtlasPositionGradientContainerType::ConstIterator  it = gradient->Begin(); 
          it != gradient->End(); ++it, ++data )
      {
      for ( int i = 0; i < 3; i++ )
        {
        *( data + i * numberOfNodes ) = it.Value()[ i ];  
        } // End loop over x,y,z coordinates
        
      } // End loop over all points
    
    }
  
protected:
  EvaluateMeshPosition() {};
  virtual ~EvaluateMeshPosition() {};


  EvaluateMeshPosition(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl


