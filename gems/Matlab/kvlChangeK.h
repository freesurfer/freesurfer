#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlAtlasMesh.h"
#include "mex.h"

namespace kvl
{

class ChangeK : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef ChangeK         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( ChangeK, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
    // kvlChangeK( mesh, scaleFactor,label )
    mexPrintf("Right one");

              
    // Make sure input arguments are correct
    if ( ( nrhs != 3 ) || ( nlhs != 0 ) || !mxIsDouble( prhs[ 1 ] ) )
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

    double scaleFactor = *( mxGetPr( prhs[ 1 ] ) );
    const int classNumber =  static_cast< int >( *( mxGetPr( prhs[ 2 ] ) ) );

    if(scaleFactor<1){
      scaleFactor=scaleFactor-1;
    }

    else if(scaleFactor>1){
      scaleFactor=scaleFactor-1;
    }
    else{
      scaleFactor=1;
    }
  
    //std::cout << "mesh: " << mesh.GetPointer() << std::endl;
    //std::cout << "scaleFactor: " << scaleFactor << std::endl;
    

    // Get the alpha from each tetrahedra node corresponding to the given label and 
    // scale the K by the scaling factor times the probability
    for ( AtlasMesh::CellsContainer::Iterator cellIt = mesh->GetCells()->Begin();
          cellIt != mesh->GetCells()->End(); ++cellIt )
      {

	AtlasMesh::CellType*  cell = cellIt.Value();
	if( cell->GetType() == AtlasMesh::CellType::TETRAHEDRON_CELL ){
	  
	  AtlasMesh::CellType::PointIdIterator  pit = cell->PointIdsBegin();
	  AtlasAlphasType alphas0  = mesh->GetPointData()->ElementAt(*pit).m_Alphas;
	  ++pit;
	  AtlasAlphasType alphas1  =  mesh->GetPointData()->ElementAt(*pit).m_Alphas;
	  ++pit;
	  AtlasAlphasType alphas2  =  mesh->GetPointData()->ElementAt(*pit).m_Alphas;
	  ++pit;
	  AtlasAlphasType alphas3  =  mesh->GetPointData()->ElementAt(*pit).m_Alphas;
	  ++pit;
	  mesh->GetCellData()->ElementAt(cellIt->Index()).m_ReferenceVolumeTimesK *= (1+(scaleFactor/4)*(alphas0[classNumber]+alphas1[classNumber]+alphas2[classNumber]+alphas3[classNumber]));
	

	}

 
      }


    }


protected:
  ChangeK() {};
  virtual ~ChangeK() {};


  ChangeK(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl



