#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlCroppedImageReader.h"
#include "kvlAtlasMeshCollection.h"
#include "vnl/vnl_det.h"


namespace kvl
{

class ReadMeshCollection : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef ReadMeshCollection         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( ReadMeshCollection, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
    // meshCollection = kvlReadMeshCollection( meshCollectionFileName, transform, K )
  
    // Retrieve the input arguments
    if ( nrhs < 1 )
      {
      mexErrMsgTxt( "meshCollectionFileName needed" );
      }
    if ( !mxIsChar( prhs[ 0 ] ) )
      {
      mexErrMsgTxt( "meshCollectionFileName must be string" );
      }
    const std::string  meshCollectionFileName = mxArrayToString( prhs[0] );
    
    typedef CroppedImageReader::TransformType  TransformType;
    TransformType::ConstPointer  transform = 0;
    float  K = -1.0f;
    
    if ( nrhs > 1 ) 
      {
      if ( !mxIsInt64( prhs[ 1 ] ) )
        {
        mexErrMsgTxt( "transform must be in int64 format" );
        }
        
      const int transformHandle = *( static_cast< int* >( mxGetData( prhs[ 1 ] ) ) );
      itk::Object::ConstPointer object = kvl::MatlabObjectArray::GetInstance()->GetObject( transformHandle );
 
      // if ( typeid( *object ) != typeid( TransformType ) )
      if ( strcmp(typeid( *object ).name(), typeid( TransformType ).name()) )  // Eugenio: MAC compatibility
        {
        mexErrMsgTxt( "transform doesn't refer to the correct ITK object type" );
        }
        
      transform = static_cast< const TransformType* >( object.GetPointer() );
      }
    
    if ( nrhs > 2 )
      {
      K = static_cast< float >( *( mxGetPr( prhs[ 2 ] ) ) );
      }
    
    //std::cout << "meshCollectionFileName: " << meshCollectionFileName << std::endl;
    //std::cout << "transform: " << transform.GetPointer() << std::endl;
    //std::cout << "K: " << K << std::endl;

    
    
    
    // Read the mesh collection
    kvl::AtlasMeshCollection::Pointer  meshCollection = kvl::AtlasMeshCollection::New();

    if ( !meshCollection->Read( meshCollectionFileName.c_str() ) )
      {
      itkExceptionMacro( "Couldn't read mesh collection from file " << meshCollectionFileName );
      }

    // Change K if user has specified a value
    if ( K > 0 )
      {
      //std::cout << "Setting K of mesh collection to: " << K << std::endl;
      meshCollection->SetK( K );
      }


    // Apply the correct transform
    if ( transform )
      {
      //std::cout << "Applying transform: " << std::endl;
      meshCollection->Transform( -1, transform );
      for ( unsigned int i = 0; i < meshCollection->GetNumberOfMeshes(); i++ )
        {
        meshCollection->Transform( i, transform );
        }

      const float  determinant = vnl_det( transform->GetMatrix().GetVnlMatrix() );
      if ( determinant < 0 )
        {
        //std::cout << "Careful here: the applied transformation will turn positive tetrahedra into negative ones." << std::endl;
        //std::cout << transform->GetMatrix().GetVnlMatrix() << std::endl;
        //std::cout << " determinant: " << determinant << std::endl;
        //std::cout << "Starting to swap the point assignments of each tetrahedron..." << std::endl;

        for ( kvl::AtlasMesh::CellsContainer::Iterator  cellIt = meshCollection->GetCells()->Begin();
              cellIt != meshCollection->GetCells()->End(); ++cellIt )
          {
          kvl::AtlasMesh::CellType*  cell = cellIt.Value();

          if( cell->GetType() != kvl::AtlasMesh::CellType::TETRAHEDRON_CELL )
            {
            continue;
            }

          // Swap points assigned to first two vertices. This will readily turn negative tetrahedra
          //into positives ones.
          kvl::AtlasMesh::CellType::PointIdIterator  pit = cell->PointIdsBegin();
          const kvl::AtlasMesh::PointIdentifier  p0Id = *pit;
          ++pit;
          const kvl::AtlasMesh::PointIdentifier  p1Id = *pit;

          pit = cell->PointIdsBegin();
          *pit = p1Id;
          ++pit;
          *pit = p0Id;
          } // End loop over all tetrahedra


        //std::cout << "...done!" << std::endl;
        }
        
      } // End test if a transform is given

      
      
    // Store the meshCollection in persistent memory
    const int meshCollectionHandle = kvl::MatlabObjectArray::GetInstance()->AddObject( meshCollection );
     
    // Return the handle to Matlab
    mwSize  dims[ 1 ];
    dims[ 0 ] = 1;
    plhs[ 0 ] = mxCreateNumericArray( 1, dims, mxINT64_CLASS, mxREAL );
    *( static_cast< int* >( mxGetData( plhs[ 0 ] ) ) ) = meshCollectionHandle;

    }



protected:
  ReadMeshCollection() {};
  virtual ~ReadMeshCollection() {};


  ReadMeshCollection(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl



