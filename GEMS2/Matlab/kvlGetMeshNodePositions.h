#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlAtlasMesh.h"


namespace kvl
{

class GetMeshNodePositions : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef GetMeshNodePositions         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( GetMeshNodePositions, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
              
    // positions = kvlGetMeshNodePositions( mesh )
  
    // Make sure input arguments are correct
    if ( ( nrhs != 1 ) || !mxIsInt64( prhs[ 0 ] ) || 
         ( nlhs != 1 ) )
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
    kvl::AtlasMesh::ConstPointer mesh = static_cast< const kvl::AtlasMesh* >( object.GetPointer() );


    // Copy the alphas in the mesh nodes into a Matlab matrix
    const int  numberOfNodes = mesh->GetPoints()->Size();
    //std::cout << "numberOfNodes :" << numberOfNodes << std::endl;
    mwSize  dims[ 2 ];
    dims[ 0 ] = numberOfNodes;
    dims[ 1 ] = 3;
    plhs[ 0 ] = mxCreateNumericArray( 2, dims, mxDOUBLE_CLASS, mxREAL );
    double*  data = static_cast< double* >( mxGetData( plhs[ 0 ] ) ); 

    for ( AtlasMesh::PointsContainer::ConstIterator  it = mesh->GetPoints()->Begin(); 
          it != mesh->GetPoints()->End(); ++it, ++data )
      {
      for ( int i = 0; i < 3; i++ )
        {
        *( data + i * numberOfNodes ) = it.Value()[ i ];  
        } // End loop over x,y,z coordinates
        
      } // End loop over all points
    
    }
  
protected:
  GetMeshNodePositions() {};
  virtual ~GetMeshNodePositions() {};


  GetMeshNodePositions(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl



