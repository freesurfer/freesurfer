#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlAtlasMesh.h"


namespace kvl
{

class SetAlphasInMeshNodes : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef SetAlphasInMeshNodes         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( SetAlphasInMeshNodes, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
              
    // kvlSetAlphasInMeshNodes( mesh, alphas )
  
    // Make sure input arguments are correct
    if ( ( nrhs != 2 ) || !mxIsInt64( prhs[ 0 ] ) || 
         !mxIsSingle( prhs[ 1 ] ) || ( mxGetNumberOfDimensions( prhs[ 1 ] ) != 2 ) || 
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
    const int  numberOfLabels = mxGetDimensions( prhs[ 1 ] )[ 1 ];
    if ( mesh->GetPointData()->Size() != numberOfNodes )
      {
      mexErrMsgTxt( "Dimensions of alphas don't match the mesh properties" );
      }
    const float*  data = static_cast< float* >( mxGetData( prhs[ 1 ] ) ); 

    
    // Copy the alphas from the Matlab matrix into the mesh nodes
    for ( AtlasMesh::PointDataContainer::Iterator  it = mesh->GetPointData()->Begin(); 
          it != mesh->GetPointData()->End(); ++it, ++data )
      {
      AtlasAlphasType  alphas( numberOfLabels );
      for ( int i = 0; i < numberOfLabels; i++ )
        {
        alphas[ i ] = *( data + i * numberOfNodes );  
        } // End loop over all labels

      it.Value().m_Alphas = alphas;

      } // End loop over all point parameters
    
    }
  
protected:
  SetAlphasInMeshNodes() {};
  virtual ~SetAlphasInMeshNodes() {};


  SetAlphasInMeshNodes(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl



