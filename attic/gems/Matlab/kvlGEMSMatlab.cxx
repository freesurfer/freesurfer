#include "mex.h" 
#include "itkMGHImageIOFactory.h"
#include "kvlMatlabRunnerArray.h"


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

   // Make sure first input argument is a string
  if ( ( nrhs == 0 ) || !mxIsChar( prhs[ 0 ] ) ) 
    {
    mexErrMsgTxt( "First argument should be string (runnerName)" );
    }

  // Retrieve the string
  const std::string runnerName = mxArrayToString( prhs[0] );

  try
    {
    // Add support for MGH file format to ITK. An alternative way to add this by default would be
    // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
    itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );

    // Try to run 
    const int nrhsToPassOn = nrhs - 1;
    const mxArray**  prhsToPassOn = prhs + 1;
    const bool  couldRun = 
       kvl::MatlabRunnerArray::GetInstance()->Run( runnerName, 
                                                   nlhs, plhs,
                                                   nrhsToPassOn, prhsToPassOn );
    if ( !couldRun )
      {
      std::ostringstream  errorStream;
      errorStream << "Coulnd't find anything to do " << runnerName;
      mexErrMsgTxt( errorStream.str().c_str() );
      }

    }
  catch ( std::exception& e )
    {
    mexErrMsgTxt( e.what() );
    return;
    }
}

