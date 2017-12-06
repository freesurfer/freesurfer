#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlAtlasMeshDeformationLevenbergMarquardtOptimizer.h"


namespace kvl
{

class GetLevenbergMarquardtOptimizer : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef GetLevenbergMarquardtOptimizer         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( GetLevenbergMarquardtOptimizer, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
              
    // optimizer = kvlGetLevenbergMarquardtOptimizer( mesh, image, transform )
  
    // Make sure input arguments are correct
    if ( ( nrhs < 3 ) || 
         !mxIsInt64( prhs[ 0 ] ) || 
         !mxIsInt64( prhs[ 1 ] ) || 
         !mxIsInt64( prhs[ 2 ] ) )
      {
      mexErrMsgTxt( "Incorrect arguments" );
      }
      
    // Retrieve input mesh
    const int meshHandle = *( static_cast< int* >( mxGetData( prhs[ 0 ] ) ) );
    itk::Object::ConstPointer object = kvl::MatlabObjectArray::GetInstance()->GetObject( meshHandle );
    // if ( typeid( *object ) != typeid( kvl::AtlasMesh ) )
    if ( strcmp(typeid( *object ).name(), typeid( kvl::AtlasMesh ).name()) )  // Eugenio: MAC compatibility
      {
      mexErrMsgTxt( "mesh doesn't refer to the correct ITK object type" );
      }
    kvl::AtlasMesh::ConstPointer constMesh = static_cast< const kvl::AtlasMesh* >( object.GetPointer() );
    kvl::AtlasMesh::Pointer mesh = const_cast< kvl::AtlasMesh* >( constMesh.GetPointer() );

    // Retrieve input image(s)
    typedef AtlasMeshDeformationLevenbergMarquardtOptimizer::ImageType  ImageType;

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
	 // if ( typeid(*(object)  ) != typeid( ImageType ) )
	 if ( strcmp(typeid( *object ).name(), typeid( ImageType ).name()) )  // Eugenio: MAC compatibility
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
    // if ( typeid( *object ) != typeid( TransformType ) )
    if ( strcmp(typeid( *object ).name(), typeid( TransformType ).name()) )  // Eugenio: MAC compatibility
      {
      mexErrMsgTxt( "transform doesn't refer to the correct ITK object type" );
      }
    TransformType::ConstPointer  constTransform = static_cast< const TransformType* >( object.GetPointer() );
    TransformType::Pointer  transform = const_cast< TransformType* >( constTransform.GetPointer() );

   
    // Show what we have so far
    //std::cout << "mesh: " << mesh.GetPointer() << std::endl;
    //std::cout << "image: " << image.GetPointer() << std::endl;
    //std::cout << "transform: " << transform.GetPointer() << std::endl;
    

    // Set up the optimizer accordingly
    AtlasMeshDeformationLevenbergMarquardtOptimizer::Pointer  optimizer = AtlasMeshDeformationLevenbergMarquardtOptimizer::New();
    optimizer->SetMeshToImageTransform( transform );
    optimizer->SetImages( images );
    optimizer->SetMesh( mesh );

    // Store the optimizer in persistent memory
    const int optimizerHandle = kvl::MatlabObjectArray::GetInstance()->AddObject( optimizer );
     
    // Return the handle to Matlab
    mwSize  dims[ 1 ];
    dims[ 0 ] = 1;
    plhs[ 0 ] = mxCreateNumericArray( 1, dims, mxINT64_CLASS, mxREAL );
    *( static_cast< int* >( mxGetData( plhs[ 0 ] ) ) ) = optimizerHandle;
    
    }
  
protected:
  GetLevenbergMarquardtOptimizer() {};
  virtual ~GetLevenbergMarquardtOptimizer() {};


  GetLevenbergMarquardtOptimizer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl


