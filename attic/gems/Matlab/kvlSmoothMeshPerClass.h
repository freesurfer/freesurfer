#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlAtlasMeshCollection.h"
#include "kvlAtlasMeshSmoother.h"


namespace kvl
{

class SmoothMesh : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef SmoothMesh         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( SmoothMesh, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
    // kvlSmoothMesh( mesh, sigma, classesToSmooth )
              
    // Make sure input arguments are correct
    if ( ( nrhs != 3 ) || !mxIsDouble( prhs[ 1 ] ) )
      {
      mexErrMsgTxt( "Incorrect arguments" );
      }
      
    // Retrieve input arguments
    const int meshHandle = *( static_cast< int* >( mxGetData( prhs[ 0 ] ) ) );
    itk::Object::ConstPointer object = kvl::MatlabObjectArray::GetInstance()->GetObject( meshHandle );
    // if ( typeid( *object ) != typeid( kvl::AtlasMesh ) )
    if ( strcmp(typeid( *object ).name(), typeid( kvl::AtlasMesh ).name()) )  // Eugenio: MAC compatibility
      {
      mexErrMsgTxt( "Not an atlas mesh object" );
      }
    kvl::AtlasMesh::ConstPointer  constMesh
            = static_cast< const kvl::AtlasMesh* >( object.GetPointer() );
    kvl::AtlasMesh::Pointer  mesh
            = const_cast< kvl::AtlasMesh* >( constMesh.GetPointer() );

    const double sigma = *( mxGetPr( prhs[ 1 ] ) );

    std::vector<int> classesToSmooth;
    
    const int  numberOfClasses = mxGetDimensions( prhs[ 2 ] )[ 0 ];

    std::cout << "Number of Classes: " << numberOfClasses << std::endl;

    for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
	classesToSmooth.push_back((mxGetPr( prhs[ 2 ] ))[classNumber]);
      }
    //std::cout << "mesh: " << mesh.GetPointer() << std::endl;
    //std::cout << "sigma: " << sigma << std::endl;

    // Construct a tempory mesh collection
    kvl::AtlasMeshCollection::Pointer  collection = kvl::AtlasMeshCollection::New();
    collection->GenerateFromSingleMesh( mesh, 1, 1000.0f );

    std::cout<<"Classes we are about to smooth are: "<<std::endl;

    for(int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
	{
	std::cout<<"classesToSmooth["<<classNumber<<"]: "<<classesToSmooth[classNumber]<<std::endl;
	}
 
    // Smooth
    kvl::AtlasMeshSmoother::Pointer  smoother = kvl::AtlasMeshSmoother::New();
    smoother->SetMeshCollection( collection );
    smoother->SetSigma( sigma );
    smoother->SetClassesToSmooth( classesToSmooth );
    kvl::AtlasMeshCollection::ConstPointer smoothedMeshCollection = smoother->GetSmoothedMeshCollection().GetPointer();

    // Set the alphas of the mesh to the smoothed version
    mesh->SetPointData( const_cast< AtlasMesh::PointDataContainer* >( smoothedMeshCollection->GetPointParameters() ) );
    }


protected:
  SmoothMesh() {};
  virtual ~SmoothMesh() {};


  SmoothMesh(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl



