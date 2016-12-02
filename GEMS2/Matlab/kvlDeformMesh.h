#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlAtlasMesh.h"
#include "kvlAtlasMeshSegmenter.h"
#include "kvlCroppedImageReader.h"


namespace kvl
{

  
  
class DeformMeshHelper : public AtlasMeshSegmenter
{
  
public:  
  
   /** Smart pointer typedef support. */
  typedef DeformMeshHelper         Self;
  typedef AtlasMeshSegmenter  Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( DeformMeshHelper, itk::Object );

  void SetMeans( std::vector< vnl_vector<float> >& means )
    {
    const_cast< std::vector< vnl_vector<float> >& >( this->GetMeans() ) = means;
    }

  void SetPrecisions( std::vector< vnl_matrix<float> >& precisions )
    {
    const_cast< std::vector< vnl_matrix<float> >& >( this->GetPrecisions() ) = precisions;
    }
      
  bool Go()
    {
    return Superclass::UpdatePosition( false );
    }
  
protected:
  DeformMeshHelper() {};
  virtual ~DeformMeshHelper() {};

  DeformMeshHelper(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};
 
  
  
  
class DeformMesh : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef DeformMesh         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( DeformMesh, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
              
    // converged = kvlDeformMesh( mesh, image, transform, means, variances, maximumNumberOfIterations, maximalDeformationStopCriterion )
  
    // Make sure input arguments are correct
    if ( ( nrhs < 7 ) || 
         !mxIsInt64( prhs[ 0 ] ) || 
         !mxIsInt64( prhs[ 1 ] ) || 
         !mxIsInt64( prhs[ 2 ] ) || 
         !mxIsDouble( prhs[ 3 ] ) ||
         !mxIsDouble( prhs[ 4 ] ) ||
         !mxIsDouble( prhs[ 5 ] ) ||
         !mxIsDouble( prhs[ 6 ] ) )
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
    typedef itk::Image< unsigned short, 3 >  ImageType;
   

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
    TransformType::ConstPointer  transform = static_cast< const TransformType* >( object.GetPointer() );
   
    
    
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
    
      
      
    // Retrieve maximumNumberOfIterations
    const int  maximumNumberOfIterations = static_cast< int >( *mxGetPr( prhs[ 5 ] )  );
      
    // Retrieve maximalDeformationStopCriterion
    const int  maximalDeformationStopCriterion = static_cast< int >( *mxGetPr( prhs[ 6 ] )  );
   
    // Show what we have so far
    std::cout << "mesh: " << mesh.GetPointer() << std::endl;
    std::cout << "image: " << images[0].GetPointer() << std::endl;
    std::cout << "transform: " << transform.GetPointer() << std::endl;
    std::cout << "means: [";
    for ( unsigned int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
	vnl_matrix<float> precMat = precisions[classNumber];
	for ( unsigned int row = 0; row < numberOfImages; row++ )
	  {
	    for ( unsigned int col = 0; col < numberOfImages; col++ )
	      {
		mexPrintf("precisions[%d][%d][%d] = %f\n",row,col,classNumber,precMat[row][col]);
	      }
	  }
      }

    for (unsigned int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
	vnl_vector<float> miini = means[classNumber];
	for ( unsigned int nima = 0; nima < numberOfImages; nima++ )
	  {
	    mexPrintf("means[%d][%d] = %f\n",nima, classNumber, miini[nima]);
	  }
       }
    std::cout << "maximumNumberOfIterations: " << maximumNumberOfIterations << std::endl;
    std::cout << "maximalDeformationStopCriterion: " << maximalDeformationStopCriterion << std::endl;
       

    // Let the beast go
    std::cout << std::flush;
    DeformMeshHelper::Pointer  helper = DeformMeshHelper::New();
    helper->SetMeshToImageTranform( transform );
    for(unsigned int nima = 0; nima < numberOfImages; nima++)
    {
      helper->SetImage(nima, images[nima] );
    }
    helper->SetMesh( mesh );
    helper->SetMeans( means );
    helper->SetPrecisions( precisions );
    helper->SetPositionUpdatingMaximumNumberOfIterations( maximumNumberOfIterations );
    helper->SetMaximalDeformationStopCriterion( maximalDeformationStopCriterion );
    const bool converged = !( helper->Go() );
    
    // Internally the AtlasMeshSegmenter generates a copy of the mesh for its calculations,
    // and works on that, so make sure we actually get the results back
    mesh->SetPoints( const_cast< AtlasMesh::PointsContainer* >( helper->GetMesh()->GetPoints() ) );
    
    
    // Return the convergence state to Matlab
    plhs[ 0 ] = mxCreateDoubleScalar( converged );
    }
  
protected:
  DeformMesh() {};
  virtual ~DeformMesh() {};


  DeformMesh(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl



