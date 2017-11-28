#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlAtlasMeshCollection.h"


namespace kvl
{

class ScaleMeshCollection : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef ScaleMeshCollection         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( ScaleMeshCollection, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
    // kvlScaleMeshCollection( meshCollection, scaleFactor )
              
    // Make sure input arguments are correct
    if ( ( nrhs != 2 ) || ( nlhs != 0 ) || !mxIsDouble( prhs[ 1 ] ) )
      {
      mexErrMsgTxt( "Incorrect arguments" );
      }
      
    // Retrieve input arguments
    const int meshCollectionHandle = *( static_cast< int* >( mxGetData( prhs[ 0 ] ) ) );
    itk::Object::ConstPointer object = kvl::MatlabObjectArray::GetInstance()->GetObject( meshCollectionHandle );
    // if ( typeid( *object ) != typeid( kvl::AtlasMeshCollection ) )
    if ( strcmp(typeid( *object ).name(), typeid( kvl::AtlasMeshCollection ).name()) )  // Eugenio: MAC compatibility
      {
      mexErrMsgTxt( "Not an atlas mesh collection object" );
      }
    kvl::AtlasMeshCollection::ConstPointer  constMeshCollection
            = static_cast< const kvl::AtlasMeshCollection* >( object.GetPointer() );
    kvl::AtlasMeshCollection::Pointer  meshCollection
            = const_cast< kvl::AtlasMeshCollection* >( constMeshCollection.GetPointer() );

    const double scaleFactor = *( mxGetPr( prhs[ 1 ] ) );
  
    //std::cout << "meshCollection: " << meshCollection.GetPointer() << std::endl;
    //std::cout << "scaleFactor: " << scaleFactor << std::endl;

    // Scale the mesh collection
    typedef kvl::AtlasMeshCollection::TransformType  TransformType;
    TransformType::Pointer  transform = TransformType::New();
    transform->Scale( scaleFactor );
    meshCollection->Transform( -1, transform );
    for ( unsigned int i = 0; i < meshCollection->GetNumberOfMeshes(); i++ )
      {
      meshCollection->Transform( i, transform );
      }

    }


protected:
  ScaleMeshCollection() {};
  virtual ~ScaleMeshCollection() {};


  ScaleMeshCollection(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl



