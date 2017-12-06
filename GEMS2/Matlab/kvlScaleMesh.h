#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlAtlasMesh.h"


namespace kvl
{

class ScaleMesh : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef ScaleMesh         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( ScaleMesh, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
    // kvlScaleMesh( mesh, scaleFactor(s) )
              
    // Make sure input arguments are correct
    if ( ( nrhs != 2 ) || ( nlhs != 0 ) || !mxIsDouble( prhs[ 1 ] ) )
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
    //std::cout << "mesh: " << mesh.GetPointer() << std::endl;

    double  scaling[ 3 ];
    double*  tmp = mxGetPr( prhs[ 1 ] );
    const int  numberOfRows = *( mxGetDimensions( prhs[ 1 ] ) );
    const int  numberOfColumns = *( mxGetDimensions( prhs[ 1 ] ) + 1 );
    std::cout << "numberOfRows: " << numberOfRows << std::endl;
    std::cout << "numberOfColumns: " << numberOfColumns << std::endl;
    if ( ( numberOfRows * numberOfColumns ) == 1 )
      {
      for ( int i = 0; i < 3; i++ )
        {
        scaling[ i ] = *tmp;
        }
        
      }
    else if ( ( numberOfRows * numberOfColumns ) == 3 )
      {
      for ( int i = 0; i < 3; i++, tmp++ )
        {
        scaling[ i ] = *tmp;
        }
      }
    else
      {
      mexErrMsgTxt( "Scale factor(s) must be 1- or 3-dimensional" );
      }


    // Transform the mesh points
    for ( AtlasMesh::PointsContainer::Iterator it = mesh->GetPoints()->Begin();
          it != mesh->GetPoints()->End(); ++it )
      {
      for ( int i = 0; i < 3; i++ )
        {
        it.Value()[i] *= scaling[ i ];
        }
      }


    // Also the reference position of the mesh has changed. Note, however, that we don't
    // actually have access to the original reference position, only some sufficient
    // statistics calculated from it. In particular, we have only the three first columns
    // of the matrix Z = inv( X ) = inv( [ p0 p1 p2 p3; 1 1 1 1 ] ). The transformed
    // position Xtrans is given by
    //
    //     Xtrans = diag( scaling[0] scaling[ 1 ] scaling[ 2 ] 1 ) * X + [ translation[ 0 ] translation[ 1 ] translation[ 2 ] 1 ]'
    //
    // but since we'll also end up calculating the upper 3x3 matrix of Ytrans * Ztrans
    // (with Y the equivalent of X but in the deformed mesh)
    // to evaluate the deformation penatly, we can safely drop the translation part.
    // In short, we will calculate the first 3 columns of the matrix
    //
    //    Ztrans = inv( diag( scaling[0] scaling[ 1 ] scaling[ 2 ] 1 ) * X )
    //
    //           = Z * diag( 1/scaling[0] 1/scaling[1] 1/scaling[2] 1 )
    //
    // which is given by multiplying each column i of Z with a factor 1/scaling[i]
    //
    for ( AtlasMesh::CellDataContainer::Iterator  it = mesh->GetCellData()->Begin();
          it != mesh->GetCellData()->End(); ++it )
      {
      it.Value().m_ReferenceVolumeTimesK *= ( scaling[ 0 ] * scaling[ 1 ] * scaling[ 2 ] );

      it.Value().m_Z11 /= scaling[ 0 ];
      it.Value().m_Z21 /= scaling[ 0 ];
      it.Value().m_Z31 /= scaling[ 0 ];
      it.Value().m_Z41 /= scaling[ 0 ];

      it.Value().m_Z12 /= scaling[ 1 ];
      it.Value().m_Z22 /= scaling[ 1 ];
      it.Value().m_Z32 /= scaling[ 1 ];
      it.Value().m_Z42 /= scaling[ 1 ];

      it.Value().m_Z13 /= scaling[ 2 ];
      it.Value().m_Z23 /= scaling[ 2 ];
      it.Value().m_Z33 /= scaling[ 2 ];
      it.Value().m_Z43 /= scaling[ 2 ];
      }


    }


protected:
  ScaleMesh() {};
  virtual ~ScaleMesh() {};


  ScaleMesh(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl



