#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"


namespace kvl
{

  
class CreateTransform : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef CreateTransform         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( CreateTransform, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
              
    // transform = kvlCreateTransform( transformMatrix )
  
    // Make sure input arguments are correct
    if ( ( nrhs != 1 ) || 
         !mxIsDouble( prhs[ 0 ] ) || //Single/>double
         ( mxGetNumberOfDimensions( prhs[ 0 ] ) != 2 ) ||
         ( mxGetDimensions( prhs[ 0 ] )[ 0 ] != 4 ) ||
         ( mxGetDimensions( prhs[ 0 ] )[ 1 ] != 4 ) ||
         ( nlhs != 1 ) )
      {
      mexErrMsgTxt( "Incorrect arguments" );
      }
      
      
    // Create the ITK transform object and fill in its elements
    typedef itk::AffineTransform< double, 3 >  TransformType; //double here!
    TransformType::Pointer  transform = TransformType::New();
    const double*  data = static_cast< const double* >( mxGetData( prhs[ 0 ] ) ); //float->double
    TransformType::ParametersType  parameters( 12 );
    for ( unsigned int row = 0; row < 3; row++ )
      {
      for ( unsigned int col = 0; col < 3; col++ )
        {
        parameters[ row * 3 + col ] = data[ col * 4 + row ];
        }
      parameters[ 9 + row ] = data[ 12 + row ];
      }       
    transform->SetParameters( parameters );  
      
    
    // Store the transform in persistent memory
    const int transformHandle = kvl::MatlabObjectArray::GetInstance()->AddObject( transform );
    
    // Return the handle to Matlab
    mwSize  dims[ 1 ];
    dims[ 0 ] = 1;
    plhs[ 0 ] = mxCreateNumericArray( 1, dims, mxINT64_CLASS, mxREAL );
    *( static_cast< int* >( mxGetData( plhs[ 0 ] ) ) ) = transformHandle;
      
    }
  
protected:
  CreateTransform() {};
  virtual ~CreateTransform() {};


  CreateTransform(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:
  


};

} // end namespace kvl


