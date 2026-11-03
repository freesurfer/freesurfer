#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlImageConverter.h"
#include "itkImage.h"
#include "itkImageRegionConstIterator.h"


namespace kvl
{

  
class GetTransformMatrix : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef GetTransformMatrix         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( GetTransformMatrix, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
              
    // transformMatrix = kvlGetTransformMatrix( transform )
  
    // Make sure input arguments are correct
    if ( ( nrhs != 1 ) || !mxIsInt64( prhs[ 0 ] ) || 
         ( nlhs != 1 ) )
      {
      mexErrMsgTxt( "Incorrect arguments" );
      }
      
    // Retrieve the tranform
    typedef itk::AffineTransform< double, 3 >  TransformType; //double here!
    const int transformHandle = *( static_cast< int* >( mxGetData( prhs[ 0 ] ) ) );
    itk::Object::ConstPointer object = kvl::MatlabObjectArray::GetInstance()->GetObject( transformHandle );
    // if ( typeid( *object ) != typeid( TransformType ) )
    if ( strcmp(typeid( *object ).name(), typeid(  TransformType ).name()) )  // Eugenio: MAC compatibility
      {
      mexErrMsgTxt( "transform doesn't refer to the correct ITK object type" );
      }
    TransformType::ConstPointer  transform = static_cast< const TransformType* >( object.GetPointer() );

    
    // Create a Matlab matrix and fill in
    mwSize  dims[ 2 ];
    dims[ 0 ] = 4;
    dims[ 1 ] = 4;
    plhs[ 0 ] = mxCreateNumericArray( 2, dims, mxDOUBLE_CLASS, mxREAL );
    double*  data = static_cast< double* >( mxGetData( plhs[ 0 ] ) ); 
    TransformType::ParametersType  parameters = transform->GetParameters();

    for ( unsigned int row = 0; row < 3; row++ )
      {
      for ( unsigned int col = 0; col < 3; col++ )
        {
        data[ col * 4 + row ] = parameters[ row * 3 + col ];
        }
      data[ 12 + row ] = parameters[ 9 + row ];
      }
    for ( unsigned int col = 0; col < 3; col++ )
      {
      data[ col * 4 + 3 ] = 0.0f;
      }
    data[ 15 ] = 1.0f;
    
    
    }
  
protected:
  GetTransformMatrix() {};
  virtual ~GetTransformMatrix() {};


  GetTransformMatrix(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:
  


};

} // end namespace kvl


