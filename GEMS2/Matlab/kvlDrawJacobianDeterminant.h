#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "itkImage.h"
#include "itkImageRegionConstIterator.h"
#include "kvlAtlasMeshJacobianDeterminantDrawer.h"



namespace kvl
{

class DrawJacobianDeterminant : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef DrawJacobianDeterminant         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( DrawJacobianDeterminant, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
              
    // jacobianDeterminantBuffer = kvlDrawJacobianDeterminant( mesh, DIM )
  
    // Make sure input arguments are correct
    if ( ( nrhs < 2 ) || !mxIsInt64( prhs[ 0 ] ) || 
         ( nlhs != 1 ) )
      {
      mexErrMsgTxt( "Incorrect arguments" );
      }
      
    // Some typedefs
    typedef AtlasMeshJacobianDeterminantDrawer::ImageType  JacobianDeterminantImageType;
    typedef JacobianDeterminantImageType::SizeType  SizeType;
   
    // Retrieve input arguments
    const int meshHandle = *( static_cast< int* >( mxGetData( prhs[ 0 ] ) ) );
    itk::Object::ConstPointer object = kvl::MatlabObjectArray::GetInstance()->GetObject( meshHandle );
    if ( typeid( *object ) != typeid( kvl::AtlasMesh ) )
      {
      mexErrMsgTxt( "mesh doesn't refer to the correct ITK object type" );
      }
    kvl::AtlasMesh::ConstPointer mesh = static_cast< const kvl::AtlasMesh* >( object.GetPointer() );
    double*  tmp = mxGetPr( prhs[ 1 ] );
    SizeType  imageSize;
    for ( int i = 0; i < 3; i++, tmp++ )
      {
      imageSize[ i ] = static_cast< int >( *tmp );
      }
      
    //std::cout << "mesh: " << mesh.GetPointer() << std::endl;
    //std::cout << "imageSize: " << imageSize << std::endl;
    

    // Rasterize
    kvl::AtlasMeshJacobianDeterminantDrawer::Pointer  drawer = kvl::AtlasMeshJacobianDeterminantDrawer::New();
    drawer->SetRegions( imageSize );
    drawer->Rasterize( mesh );

    
    // Finally, copy the buffer contents into a Matlab matrix
    mwSize  dims[ 3 ];
    for ( int i = 0; i < 3; i++ )
      {
      dims[ i ] = imageSize[ i ];
      }
    //plhs[ 0 ] = mxCreateNumericArray( 3, dims, mxSINGLE_CLASS, mxREAL );
    //float*  data = static_cast< float* >( mxGetData( plhs[ 0 ] ) ); 
    plhs[ 0 ] = mxCreateNumericArray( 3, dims, mxSINGLE_CLASS, mxREAL );
    float*  data = static_cast< float* >( mxGetData( plhs[ 0 ] ) ); 

    itk::ImageRegionConstIterator< JacobianDeterminantImageType >  
                      it( drawer->GetImage(),
                          drawer->GetImage()->GetBufferedRegion() );
    for ( ;!it.IsAtEnd(); ++it, ++data )
      {
      *data = it.Value();
      }

    }
  
protected:
  DrawJacobianDeterminant() {};
  virtual ~DrawJacobianDeterminant() {};


  DrawJacobianDeterminant(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl



