#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlAtlasMeshCollection.h"
#include "vnl/vnl_det.h"


namespace kvl
{

class GetMesh : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef GetMesh         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( GetMesh, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
    // mesh = kvlGetMesh( meshCollection, meshNumber )
  
    if ( ( nrhs < 1 ) || ( nlhs != 1 ) )
      {
      mexErrMsgTxt( "Wrong number of input/output arguments" );
      }
              
    // Retrieve input arguments
    const int meshCollectionHandle = *( static_cast< int* >( mxGetData( prhs[ 0 ] ) ) );
    itk::Object::ConstPointer object = kvl::MatlabObjectArray::GetInstance()->GetObject( meshCollectionHandle );
    // if ( typeid( *object ) != typeid( kvl::AtlasMeshCollection ) )
    if ( strcmp(typeid( *object ).name(), typeid( kvl::AtlasMeshCollection ).name()) )  // Eugenio: MAC compatibility
      {
      mexErrMsgTxt( "Not an atlas mesh collection object" );
      }
    kvl::AtlasMeshCollection::ConstPointer  meshCollection
            = static_cast< const kvl::AtlasMeshCollection* >( object.GetPointer() );
  
    int  meshNumber = -1;
    if ( nrhs > 1 )
      {
      meshNumber = static_cast< int >( *( mxGetPr( prhs[ 1 ] ) ) );
      }
    
    //std::cout << "meshCollection: " << meshCollection.GetPointer() << std::endl;
    //std::cout << "meshNumber: " << meshNumber << std::endl;


    // Now get the correct mesh
    kvl::AtlasMesh::ConstPointer  mesh = 0;
    if ( meshNumber == -1 )
      {
      //std::cout << "Getting reference mesh" << std::endl;
      mesh = meshCollection->GetReferenceMesh();
      }
    else
      {
      //std::cout << "Getting mesh number: " << meshNumber << std::endl;
      mesh = meshCollection->GetMesh( meshNumber );
      }
      
      
    // Store the mesh in persistent memory
    const int meshHandle = kvl::MatlabObjectArray::GetInstance()->AddObject( mesh );
     
    // Return the handle to Matlab
    mwSize  dims[ 1 ];
    dims[ 0 ] = 1;
    plhs[ 0 ] = mxCreateNumericArray( 1, dims, mxINT64_CLASS, mxREAL );
    *( static_cast< int* >( mxGetData( plhs[ 0 ] ) ) ) = meshHandle;

    }



protected:
  GetMesh() {};
  virtual ~GetMesh() {};


  GetMesh(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl



