#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"


namespace kvl
{

  
class Clear: public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef Clear         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( Clear, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
    // kvlClear( handle )
    
    // Clear all if no handle is given
    if ( nrhs == 0 )
      {
      kvl::MatlabObjectArray::GetInstance()->Clear();
      return;
      }


    if ( !mxIsInt64( prhs[ 0 ] ) )
      {
      mexErrMsgTxt( "Incorrect arguments" );
      }
      
    // Retrieve input handle
    const int  handle = *( static_cast< int* >( mxGetData( prhs[ 0 ] ) ) );
    //std::cout << "handle: " << handle << std::endl;
    kvl::MatlabObjectArray::GetInstance()->RemoveObject( handle );
    }
  
protected:
  Clear() {};
  virtual ~Clear() {};


  Clear(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl



