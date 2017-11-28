#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlAtlasMesh.h"


namespace kvl
{

class GetAlphasInMeshNodes : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef GetAlphasInMeshNodes         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( GetAlphasInMeshNodes, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
              
    // alphas = kvlGetAlphasInMeshNodes( mesh )
  
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
    const int  numberOfNodes = mesh->GetPointData()->Size();
    const int  numberOfLabels = mesh->GetPointData()->Begin().Value().m_Alphas.Size();
    //std::cout << "numberOfNodes :" << numberOfNodes << std::endl;
    //std::cout << "numberOfLabels:" << numberOfLabels << std::endl;
    mwSize  dims[ 2 ];
    dims[ 0 ] = numberOfNodes;
    dims[ 1 ] = numberOfLabels;
    plhs[ 0 ] = mxCreateNumericArray( 2, dims, mxSINGLE_CLASS, mxREAL );
    float*  data = static_cast< float* >( mxGetData( plhs[ 0 ] ) ); 

    for ( AtlasMesh::PointDataContainer::ConstIterator  it = mesh->GetPointData()->Begin(); 
          it != mesh->GetPointData()->End(); ++it, ++data )
      {
      for ( int i = 0; i < numberOfLabels; i++ )
        {
        *( data + i * numberOfNodes ) = it.Value().m_Alphas[ i ];  
        } // End loop over all labels
        
      } // End loop over all point parameters
    
    }
  
protected:
  GetAlphasInMeshNodes() {};
  virtual ~GetAlphasInMeshNodes() {};


  GetAlphasInMeshNodes(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl



