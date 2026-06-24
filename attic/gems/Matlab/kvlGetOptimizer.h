#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlAtlasMeshDeformationFixedStepGradientDescentOptimizer.h"
#include "kvlAtlasMeshDeformationGradientDescentOptimizer.h"
#include "kvlAtlasMeshDeformationConjugateGradientOptimizer.h"
#include "kvlAtlasMeshDeformationLBFGSOptimizer.h"


namespace kvl
{

class GetOptimizer : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef GetOptimizer          Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( GetOptimizer, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
              
    // Make sure input arguments are correct
    const  std::string  usageString = "Usage: optimizer = kvlGetOptimizer( typeName, mesh, calculator )\n where typeName = {'FixedStepGradientDescent','GradientDescent','ConjugateGradient','L-BFGS'}";
    if ( ( nrhs < 3 ) || 
         !mxIsChar( prhs[ 0 ] ) || 
         !mxIsInt64( prhs[ 1 ] ) ||
         !mxIsInt64( prhs[ 2 ] ) ) 
      {
      mexErrMsgTxt( usageString.c_str() ); 
      }

    // Retrieve the mesh
    const int meshHandle = *( static_cast< int* >( mxGetData( prhs[ 1 ] ) ) );
    itk::Object::ConstPointer  object = kvl::MatlabObjectArray::GetInstance()->GetObject( meshHandle );
    kvl::AtlasMesh::ConstPointer  constMesh
             = dynamic_cast< const kvl::AtlasMesh* >( object.GetPointer() );
    if ( !constMesh.GetPointer() )
      {
      std::cout << "typeid: " << typeid( *object ).name() << std::endl;
      mexErrMsgTxt( "mesh doesn't refer to the correct ITK object type" );
      }  
    kvl::AtlasMesh::Pointer  mesh = const_cast< kvl::AtlasMesh* >( constMesh.GetPointer() );
      
                   
    // Retrieve the calculator
    const int calculatorHandle = *( static_cast< int* >( mxGetData( prhs[ 2 ] ) ) );
    object = kvl::MatlabObjectArray::GetInstance()->GetObject( calculatorHandle );
    kvl::AtlasMeshPositionCostAndGradientCalculator::ConstPointer  constCalculator
             = dynamic_cast< const kvl::AtlasMeshPositionCostAndGradientCalculator* >( object.GetPointer() );
    if ( !constCalculator.GetPointer() )
      {
      std::cout << "typeid: " << typeid( *object ).name() << std::endl;
      mexErrMsgTxt( "calculator doesn't refer to the correct ITK object type" );
      }  
    kvl::AtlasMeshPositionCostAndGradientCalculator::Pointer  calculator 
                   = const_cast< kvl::AtlasMeshPositionCostAndGradientCalculator* >( constCalculator.GetPointer() );

                   
    // Construct the correct type of optimizer
    AtlasMeshDeformationOptimizer::Pointer  optimizer = 0; 
    const std::string  typeName = mxArrayToString( prhs[ 0 ] );
    switch( typeName[ 0 ] ) 
      {
      case 'F': 
        {
        std::cout << "FixedStepGradientDescent" << std::endl;
        AtlasMeshDeformationFixedStepGradientDescentOptimizer::Pointer  myOptimizer
                          = AtlasMeshDeformationFixedStepGradientDescentOptimizer::New();
        myOptimizer->SetStepSize( 1.0 );
        optimizer = myOptimizer;
        break;
        } 
      case 'G': 
        {
        std::cout << "GradientDescent" << std::endl;
        AtlasMeshDeformationGradientDescentOptimizer::Pointer  myOptimizer
                          = AtlasMeshDeformationGradientDescentOptimizer::New();
        optimizer = myOptimizer;
        break;
        } 
      case 'C': 
        {
        std::cout << "ConjugateGradient" << std::endl;
        AtlasMeshDeformationConjugateGradientOptimizer::Pointer  myOptimizer 
                          = AtlasMeshDeformationConjugateGradientOptimizer::New();
        optimizer = myOptimizer;
        break;
        } 
      case 'L': 
        {
        std::cout << "L-BFGS" << std::endl;
        AtlasMeshDeformationLBFGSOptimizer::Pointer  myOptimizer 
                          = AtlasMeshDeformationLBFGSOptimizer::New();
        optimizer = myOptimizer;
        break;
        } 
      default:
        {
        mexErrMsgTxt( "typeName not understood" );
        break;
        }
      }
    
    
    // Parse additional options. Format is always [ 'someString', double ]
    for ( int  argumentNumber = 0; argumentNumber < (nrhs-3)/2; argumentNumber++ )
      {
      if ( !mxIsChar( prhs[ 3 + 2*argumentNumber ] ) || !mxIsDouble( prhs[ 3 + 2*argumentNumber + 1 ] ) )
        {
        mexErrMsgTxt( "Error parsing options" ); 
        }
        
      const std::string  optionName = mxArrayToString( prhs[ 3 + 2*argumentNumber ] );
      const double  optionValue = *( mxGetPr( prhs[ 3 + 2*argumentNumber + 1 ] ) );
      //std::cout << "optionName: " << optionName << std::endl;
      //std::cout << "optionValue: " << optionValue << std::endl;
    
      switch( optionName[ 0 ] ) 
        {
        case 'V': 
          {
          std::cout << "Verbose: " << optionValue << std::endl;
          if ( optionValue )
            {
            optimizer->SetVerbose( true );
            }

          break;
          } 
        case 'M': 
          {
          if ( optionName.substr( 0, 8 ) == "MaximalD" )
            {
            std::cout << "MaximalDeformationStopCriterion: " << optionValue << std::endl;
            optimizer->SetMaximalDeformationStopCriterion( optionValue );
            }
          else if ( optionName.substr( 0, 8 ) == "MaximumN" ) 
            {
            const int  maximumNumberOfIterations = static_cast< int >( optionValue );
            std::cout << "MaximumNumberOfIterations: " << maximumNumberOfIterations << std::endl;
            optimizer->SetMaximumNumberOfIterations( maximumNumberOfIterations );
            }
          else
            {
            std::ostringstream  errorStream;
            errorStream << "optionName: " << optionName << " not understood"; 
            mexErrMsgTxt( errorStream.str().c_str() );
            } // End figuring out which "M" you mean  
          
          break;
          }
        case 'L': 
          {
          std::cout << "LineSearchMaximalDeformationIntervalStopCriterion: " << optionValue << std::endl;
          optimizer->SetLineSearchMaximalDeformationIntervalStopCriterion( optionValue );
          break;
          }
        case 'B': 
          {
          AtlasMeshDeformationLBFGSOptimizer::Pointer  myOptimizer 
                 = dynamic_cast< AtlasMeshDeformationLBFGSOptimizer* >( optimizer.GetPointer() );
          if ( myOptimizer )
            {
            const int  maximumMemoryLength = static_cast< int >( optionValue );  
            std::cout << "BFGS-MaximumMemoryLength: " << maximumMemoryLength << std::endl;
            myOptimizer->SetMaximumMemoryLength( maximumMemoryLength );
            }
          else
            {
            std::cout << "BFGS-MaximumMemoryLength only applies to BFGS optimizer" << std::endl;  
            }  

          break;
          }
        default:
          {
          std::ostringstream  errorStream;
          errorStream << "optionName: " << optionName << " not understood"; 
          mexErrMsgTxt( errorStream.str().c_str() );
          break;
          }
        }
        
      }

    
    
    // Pass the mesh and calculator to it
    optimizer->SetMesh( mesh );
    optimizer->SetCostAndGradientCalculator( calculator );
    
                   
    // Store the optimizer in persistent memory
    const int optimizerHandle = kvl::MatlabObjectArray::GetInstance()->AddObject( optimizer );
     
    // Return the handle to Matlab
    mwSize  dims[ 1 ];
    dims[ 0 ] = 1;
    plhs[ 0 ] = mxCreateNumericArray( 1, dims, mxINT64_CLASS, mxREAL );
    *( static_cast< int* >( mxGetData( plhs[ 0 ] ) ) ) = optimizerHandle;
    
    }
  
protected:
  GetOptimizer() {};
  virtual ~GetOptimizer() {};


  GetOptimizer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl


