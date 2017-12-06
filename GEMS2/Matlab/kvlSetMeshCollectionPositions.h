#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlAtlasMesh.h"


namespace kvl
{

class SetMeshCollectionPositions : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef SetMeshCollectionPositions         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( SetMeshCollectionPositions, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
              
    // kvlSetMeshCollectionPositions( meshCollection, referencePosition, position0, position1, ... )
  
    // Make sure input arguments are correct
    if ( ( nrhs < 3 ) || !mxIsInt64( prhs[ 0 ] ) || 
         !mxIsDouble( prhs[ 1 ] ) || ( mxGetNumberOfDimensions( prhs[ 1 ] ) != 2 ) || 
         ( nlhs != 0 ) )
      {
      mexErrMsgTxt( "Incorrect arguments" );
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
    
    
#if 0    
    // Make sure we're not trying to write positions of non-existing meshes in the collection
    if ( meshCollection->GetNumberOfMeshes() < ( nrhs-2 ) )
      {
      mexErrMsgTxt( "More positions than there are meshes in the mesh collection" );
      }
    
    // Loop over all positions
    for ( int meshNumber = -1; meshNumber < ( nrhs-2 ); ++meshNumber )
      {
        
      // Get pointer to the correct ITK position container
      AtlasMeshCollection::PointsContainerType::Pointer  position = 0;
      if ( meshNumber < 0 )
        {
        position = meshCollection->GetReferencePosition();
        }
      else
        {
        position = meshCollection->GetPositions()[ meshNumber ];
        }  
        

      // Get pointer to the Matlab data
      const int  numberOfNodes = mxGetDimensions( prhs[ meshNumber+2 ] )[ 0 ];
      if ( position->Size() != numberOfNodes )
        {
        mexErrMsgTxt( "Number of nodes don't match" );
        }
      const double*  data = static_cast< double* >( mxGetData( prhs[ meshNumber+2 ] ) ); 

      
      // Copy the alphas from the Matlab matrix into the mesh nodes
      for ( AtlasMesh::PointsContainer::Iterator  it = position->Begin(); 
            it != position->End(); ++it, ++data )
        {
        AtlasMesh::PointType  point;
        for ( int i = 0; i < 3; i++ )
          {
          point[ i ] = *( data + i * numberOfNodes );  
          } // End loop over x,y,z directions

        it.Value() = point;

        } // End loop over all points
        
      } // End loop over all positions 
      
#else
    
    // Loop over all input position matrices, copy their content into the correct format, 
    // and save
    AtlasMeshCollection::PointsContainerType::Pointer  referencePosition = 0;
    std::vector< AtlasMeshCollection::PointsContainerType::Pointer >  positions;
    for ( int meshNumber = -1; meshNumber < ( nrhs-2 ); ++meshNumber )
      {
      // Get pointer to the Matlab data
      const int  numberOfNodes = mxGetDimensions( prhs[ meshNumber+2 ] )[ 0 ];
      if ( numberOfNodes != meshCollection->GetReferencePosition()->Size() )
        {
        mexErrMsgTxt( "Number of nodes don't match" );
        }
      const double*  data = static_cast< double* >( mxGetData( prhs[ meshNumber+2 ] ) ); 

      // Copy the coordinates from the Matlab matrix into a mesh node position container
      AtlasMeshCollection::PointsContainerType::Pointer  position 
                                           = AtlasMeshCollection::PointsContainerType::New();

      for ( AtlasMesh::PointsContainer::ConstIterator  it = meshCollection->GetReferencePosition()->Begin(); 
            it != meshCollection->GetReferencePosition()->End(); ++it, ++data )
        {
        AtlasMesh::PointType  point;
        for ( int i = 0; i < 3; i++ )
          {
          point[ i ] = *( data + i * numberOfNodes );  
          } // End loop over x,y,z directions

        position->InsertElement( it.Index(), point );
        } // End loop over all points
        
      // Save the mesh node position container for later use
      if ( meshNumber < 0 )
        {
        referencePosition = position;
        }
      else
        {
        positions.push_back( position );
        }  
        

        
      } // End loop over all input positions matrices
    
    
    // Now set the reference position and and mesh positions
    meshCollection->SetReferencePosition( referencePosition );
    meshCollection->SetPositions( positions );
    
    
#endif
    
    }
  
protected:
  SetMeshCollectionPositions() {};
  virtual ~SetMeshCollectionPositions() {};


  SetMeshCollectionPositions(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl



