#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlAtlasMeshCollection.h"


namespace kvl
{

class WriteMeshCollection : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef WriteMeshCollection         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( WriteMeshCollection, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
    // kvlWriteMeshCollection( meshCollection, fileName )
  
    // Check input arguments
    if ( ( nrhs != 2 ) || !mxIsInt64( prhs[ 0 ] ) || !mxIsChar( prhs[ 1 ] ) )
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
    AtlasMeshCollection::ConstPointer  meshCollection
            = static_cast< const AtlasMeshCollection* >( object.GetPointer() );
 
    // Retrieve file name
    const std::string  fileName = mxArrayToString( prhs[1] );

    // Write out
    meshCollection->Write( fileName.c_str() );
    
    }



protected:
  WriteMeshCollection() {};
  virtual ~WriteMeshCollection() {};


  WriteMeshCollection(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl



