#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlAtlasMeshCollection.h"
#include "kvlAtlasMeshSmoother.h"


namespace kvl
{

class SmoothMeshCollection : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef SmoothMeshCollection         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( SmoothMeshCollection, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
    // smoothedMeshCollection = kvlSmoothMeshCollection( meshCollection, sigma )
              
    // Make sure input arguments are correct
    if ( ( nrhs != 2 ) || ( nlhs != 1 ) || !mxIsDouble( prhs[ 1 ] ) )
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
    kvl::AtlasMeshCollection::ConstPointer  meshCollection
            = static_cast< const kvl::AtlasMeshCollection* >( object.GetPointer() );

    const double sigma = *( mxGetPr( prhs[ 1 ] ) );
  
    //std::cout << "meshCollection: " << meshCollection.GetPointer() << std::endl;
    //std::cout << "sigma: " << sigma << std::endl;

    //
    kvl::AtlasMeshSmoother::Pointer  smoother = kvl::AtlasMeshSmoother::New();
    smoother->SetMeshCollection( const_cast< kvl::AtlasMeshCollection* >( meshCollection.GetPointer() ) );
    smoother->SetSigma( sigma );
    kvl::AtlasMeshCollection::ConstPointer smoothedMeshCollection = smoother->GetSmoothedMeshCollection().GetPointer();

      
    // Store the smoothed mesh collection in persistent memory
    const int  smoothedMeshCollectionHandle = kvl::MatlabObjectArray::GetInstance()->AddObject( smoothedMeshCollection );
     
    // Return the handle to Matlab
    mwSize  dims[ 1 ];
    dims[ 0 ] = 1;
    plhs[ 0 ] = mxCreateNumericArray( 1, dims, mxINT64_CLASS, mxREAL );
    *( static_cast< int* >( mxGetData( plhs[ 0 ] ) ) ) = smoothedMeshCollectionHandle;

    }


protected:
  SmoothMeshCollection() {};
  virtual ~SmoothMeshCollection() {};


  SmoothMeshCollection(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl



