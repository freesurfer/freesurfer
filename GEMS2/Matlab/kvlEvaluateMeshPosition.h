#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlAtlasMeshPositionCostAndGradientCalculator.h"


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
              
              
    // [ cost gradient ] = kvlEvaluateMeshPosition( calculator, mesh ) 
  
    // Make sure input arguments are correct
     const  std::string  usageString = "Usage: [ cost gradient ] = kvlEvaluateMeshPosition( calculator, mesh )";
     if ( ( nrhs < 2 ) || 
         !mxIsInt64( prhs[ 0 ] ) || 
         !mxIsInt64( prhs[ 1 ] ) )
      {
      mexErrMsgTxt( usageString.c_str() );
      }
      
      
    // Retrieve calculator
    const int calculatorHandle = *( static_cast< int* >( mxGetData( prhs[ 0 ] ) ) );
    itk::Object::ConstPointer  object = kvl::MatlabObjectArray::GetInstance()->GetObject( calculatorHandle );
    // if ( typeid( *object ) != typeid( kvl::AtlasMeshPositionCostAndGradientCalculator ) )
    //   {
    //   std::cout << "typeid: " << typeid( *object ).name() << std::endl;
    //   mexErrMsgTxt( "calculator doesn't refer to the correct ITK object type" );
    //   }
    // kvl::AtlasMeshPositionCostAndGradientCalculator::ConstPointer  constCalculator 
    //                 = static_cast< const kvl::AtlasMeshPositionCostAndGradientCalculator* >( object.GetPointer() );
    kvl::AtlasMeshPositionCostAndGradientCalculator::ConstPointer  constCalculator
             = dynamic_cast< const kvl::AtlasMeshPositionCostAndGradientCalculator* >( object.GetPointer() );
    if ( !constCalculator.GetPointer() )
      {
      std::cout << "typeid: " << typeid( *object ).name() << std::endl;
      mexErrMsgTxt( "calculator doesn't refer to the correct ITK object type" );
      }  
    kvl::AtlasMeshPositionCostAndGradientCalculator::Pointer  calculator 
                   = const_cast< kvl::AtlasMeshPositionCostAndGradientCalculator* >( constCalculator.GetPointer() );

                   
    // Retrieve mesh
    const int meshHandle = *( static_cast< int* >( mxGetData( prhs[ 1 ] ) ) );
    object = kvl::MatlabObjectArray::GetInstance()->GetObject( meshHandle );
    // if ( typeid( *object ) != typeid( kvl::AtlasMesh ) )
    if ( strcmp(typeid( *object ).name(), typeid( kvl::AtlasMesh ).name()) )  // Eugenio: MAC compatibility
      {
      mexErrMsgTxt( "mesh doesn't refer to the correct ITK object type" );
      }
    kvl::AtlasMesh::ConstPointer constMesh = static_cast< const kvl::AtlasMesh* >( object.GetPointer() );
    //kvl::AtlasMesh::Pointer mesh = const_cast< kvl::AtlasMesh* >( constMesh.GetPointer() );


    // Let the beast go
    calculator->Rasterize( constMesh );

    // Retrieve the result
    const double  cost = calculator->GetMinLogLikelihoodTimesPrior();
    AtlasPositionGradientContainerType::ConstPointer  gradient = calculator->GetPositionGradient();
    
    
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


