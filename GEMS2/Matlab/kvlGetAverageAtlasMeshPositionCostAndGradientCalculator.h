#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlAverageAtlasMeshPositionCostAndGradientCalculator.h"


namespace kvl
{

class GetAverageAtlasMeshPositionCostAndGradientCalculator : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef GetAverageAtlasMeshPositionCostAndGradientCalculator         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( GetAverageAtlasMeshPositionCostAndGradientCalculator, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
              
    // calculator = kvlGetAverageAtlasMeshPositionCostAndGradientCalculator( meshCollection, K0, K1, transform )
  
    // Make sure input arguments are correct
    const  std::string  usageString = "Usage: calculator = kvlGetAverageAtlasMeshPositionCostAndGradientCalculator( meshCollection, K0, K1, transform )";
    if ( ( nrhs < 4 ) || 
         !mxIsInt64( prhs[ 0 ] ) ||
         !mxIsInt64( prhs[ 3 ] ) )
      {
      mexErrMsgTxt( usageString.c_str() ); 
      }

      
    // Retrieve mesh collection
    const int meshCollectionHandle = *( static_cast< int* >( mxGetData( prhs[ 0 ] ) ) );
    itk::Object::ConstPointer object = kvl::MatlabObjectArray::GetInstance()->GetObject( meshCollectionHandle );
    // if ( typeid( *object ) != typeid( kvl::AtlasMeshCollection ) )
    if ( strcmp(typeid( *object ).name(), typeid( kvl::AtlasMeshCollection ).name()) )  // Eugenio: MAC compatibility
      {
      mexErrMsgTxt( "Not an atlas mesh collection object" );
      }
    AtlasMeshCollection::ConstPointer  constMeshCollection
            = static_cast< const AtlasMeshCollection* >( object.GetPointer() );
    AtlasMeshCollection::Pointer  meshCollection = const_cast< kvl::AtlasMeshCollection* >( constMeshCollection.GetPointer() );
    
    
    // Retrieve K0
    const float  K0 = static_cast< float >( *( mxGetPr( prhs[ 1 ] ) ) );

    // Retrieve K1
    const float  K1 = static_cast< float >( *( mxGetPr( prhs[ 2 ] ) ) );


    // Retrieve transform 
    typedef CroppedImageReader::TransformType  TransformType;
    const int transformHandle = *( static_cast< int* >( mxGetData( prhs[ 3 ] ) ) );
    object = kvl::MatlabObjectArray::GetInstance()->GetObject( transformHandle );
    // if ( typeid( *object ) != typeid( TransformType ) )
    if ( strcmp(typeid( *object ).name(), typeid( TransformType ).name()) )  // Eugenio: MAC compatibility
      {
      mexErrMsgTxt( "transform doesn't refer to the correct ITK object type" );
      }
    TransformType::ConstPointer  constTransform = static_cast< const TransformType* >( object.GetPointer() );
      
    
    // Set up the calculator
    AverageAtlasMeshPositionCostAndGradientCalculator::Pointer  calculator 
                                         = AverageAtlasMeshPositionCostAndGradientCalculator::New(); 
    calculator->SetBoundaryCondition( AtlasMeshPositionCostAndGradientCalculator::SLIDING );
    calculator->SetMeshToImageTransform( constTransform );
    std::vector< AtlasMesh::PointsContainer::ConstPointer >  positions;
    std::vector< double >  Ks;   
    positions.push_back( meshCollection->GetReferencePosition() );
    Ks.push_back( K0 );
    for ( int  meshNumber = 0; meshNumber < meshCollection->GetPositions().size(); meshNumber++ )
      {
      positions.push_back( meshCollection->GetPositions()[ meshNumber ].GetPointer() );      
      Ks.push_back( K1 );  
      }
    calculator->SetPositionsAndKs( positions, Ks );
                             
    
    // Store the calculator in persistent memory
    const int calculatorHandle = kvl::MatlabObjectArray::GetInstance()->AddObject( calculator );
     
    // Return the handle to Matlab
    mwSize  dims[ 1 ];
    dims[ 0 ] = 1;
    plhs[ 0 ] = mxCreateNumericArray( 1, dims, mxINT64_CLASS, mxREAL );
    *( static_cast< int* >( mxGetData( plhs[ 0 ] ) ) ) = calculatorHandle;
    
    }
  
protected:
  GetAverageAtlasMeshPositionCostAndGradientCalculator() {};
  virtual ~GetAverageAtlasMeshPositionCostAndGradientCalculator() {};


  GetAverageAtlasMeshPositionCostAndGradientCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl


