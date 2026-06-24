#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlAtlasMesh.h"


namespace kvl
{

class SetMeshNodePositions : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef SetMeshNodePositions         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( SetMeshNodePositions, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
              
    // kvlSetMeshNodePositions( mesh, positions )
  
    // Make sure input arguments are correct
    if ( ( nrhs != 2 ) || !mxIsInt64( prhs[ 0 ] ) || 
         !mxIsDouble( prhs[ 1 ] ) || ( mxGetNumberOfDimensions( prhs[ 1 ] ) != 2 ) || 
         ( nlhs != 0 ) )
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

    
    // Get pointer to the Matlab data
    const int  numberOfNodes = mxGetDimensions( prhs[ 1 ] )[ 0 ];
    const int  numberOfDimensions = mxGetDimensions( prhs[ 1 ] )[ 1 ];
    if ( ( mesh->GetPoints()->Size() != numberOfNodes ) ||
         ( numberOfDimensions != 3 ) )
      {
      mexErrMsgTxt( "Dimensions of positions don't match the mesh properties" );
      }
    const double*  data = static_cast< double* >( mxGetData( prhs[ 1 ] ) ); 

    // Copy the positions from the Matlab matrix into the mesh nodes
    for ( AtlasMesh::PointsContainer::Iterator  it = mesh->GetPoints()->Begin(); 
          it != mesh->GetPoints()->End(); ++it, ++data )
      {
      for ( int i = 0; i < 3; i++ )
        {
        it.Value()[ i ] = *( data + i * numberOfNodes );  
        } // End loop over x,y,z coordinates
      } // End loop over vertices
    
    }
  
protected:
  SetMeshNodePositions() {};
  virtual ~SetMeshNodePositions() {};


  SetMeshNodePositions(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl



