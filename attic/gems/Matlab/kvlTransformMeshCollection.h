#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlCroppedImageReader.h"
#include "kvlAtlasMeshCollection.h"
#include "vnl/vnl_det.h"


namespace kvl
{

class TransformMeshCollection : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef TransformMeshCollection         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( TransformMeshCollection, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
    // kvlTransformMeshCollection( meshCollection, transform )
  
    // Check the input arguments
    if ( ( nrhs != 2 ) || !mxIsInt64( prhs[ 0 ] ) || !mxIsInt64( prhs[ 1 ] ) )
      {
      mexErrMsgTxt( "Incorrect arguments" );
      }


    // Retrieve the mesh collection
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


    // Retrieve the transform
    typedef CroppedImageReader::TransformType  TransformType;
    const int transformHandle = *( static_cast< int* >( mxGetData( prhs[ 1 ] ) ) );
    object = kvl::MatlabObjectArray::GetInstance()->GetObject( transformHandle );
    // if ( typeid( *object ) != typeid( TransformType ) )
    if ( strcmp(typeid( *object ).name(), typeid( TransformType ).name()) )  // Eugenio: MAC compatibility
      {
      mexErrMsgTxt( "transform doesn't refer to the correct ITK object type" );
      }
    TransformType::ConstPointer  transform = static_cast< const TransformType* >( object.GetPointer() );
    
    //std::cout << "meshCollection: " << meshCollection.GetPointer() << std::endl;
    //std::cout << "transform: " << transform.GetPointer() << std::endl;


    // Apply the transform
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

      } // End test if determinant < 0
      
    //std::cout << "...done!" << std::endl;

    }



protected:
  TransformMeshCollection() {};
  virtual ~TransformMeshCollection() {};


  TransformMeshCollection(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl



