#include "kvlAtlasMeshToIntensityImageGradientCalculatorCPU.h"

#include <itkMath.h>
#include "vnl/vnl_matrix_fixed.h"
#include "kvlTetrahedronInteriorConstIterator.h"



namespace kvl
{

//
//
//
AtlasMeshToIntensityImageGradientCalculatorCPU
::AtlasMeshToIntensityImageGradientCalculatorCPU()
{

  m_LikelihoodFilter = LikelihoodFilterType::New();  
  m_MinLogLikelihoodTimesPrior = 0;
  m_IgnoreDeformationPrior = false;
  m_OnlyDeformationPrior = false;
  m_PositionGradient = 0;
  m_Abort = false;

  this->SetMeshToImageTransform( TransformType::New() );
  
}


//
//
//
AtlasMeshToIntensityImageGradientCalculatorCPU
::~AtlasMeshToIntensityImageGradientCalculatorCPU()
{
}



//
//
//
void 
AtlasMeshToIntensityImageGradientCalculatorCPU
::SetImages( const std::vector< ImageType::Pointer >& images )
{
  // 
  for ( int contrastNumber = 0; contrastNumber < images.size(); contrastNumber++ )
    {
    m_LikelihoodFilter->SetInput( contrastNumber, images[ contrastNumber ] );
    }

}


//
//
//
void 
AtlasMeshToIntensityImageGradientCalculatorCPU
::SetMeans( const std::vector< vnl_vector<float> >& means )
{
  m_LikelihoodFilter->SetMeans( means );
  
}


//
//
//
void
AtlasMeshToIntensityImageGradientCalculatorCPU
::SetPrecisions( const std::vector< vnl_matrix<float> >& precisions )
{
  m_LikelihoodFilter->SetPrecisions( precisions );
}


//
//
//
void
AtlasMeshToIntensityImageGradientCalculatorCPU
::Rasterize( const AtlasMesh* mesh )
{

  // Make sure the likelihoods are up-to-date
  //m_LikelihoodFilter->SetNumberOfThreads( 1 );
  m_LikelihoodFilter->Update();

  // Initialize from a clean slate
  m_Abort = false;
  m_PositionGradient = 0;
  m_MinLogLikelihoodTimesPrior = 0;
  m_ThreadSpecificPositionGradients.clear();
  m_ThreadSpecificMinLogLikelihoodTimesPriors.clear();

  // For each thread, create an empty gradient and cost so that
  // different threads never interfere with one another
  for ( int threadNumber = 0; threadNumber < this->GetNumberOfThreads(); threadNumber++ )
    {
    // Initialize cost to zero for this thread
    m_ThreadSpecificMinLogLikelihoodTimesPriors.push_back( 0.0 );  
      
    // Create a container to hold the position gradient of this thread, and initialize to zero
    AtlasPositionGradientContainerType::Pointer  positionGradient = AtlasPositionGradientContainerType::New();
    AtlasPositionGradientType  zeroEntry( 0.0f );
    for ( AtlasMesh::PointsContainer::ConstIterator pointIt = mesh->GetPoints()->Begin();
          pointIt != mesh->GetPoints()->End(); ++pointIt )
      {
      positionGradient->InsertElement( pointIt.Index(), zeroEntry );
      }
    m_ThreadSpecificPositionGradients.push_back( positionGradient );
    
    } // End loop over threads
    
  // Now rasterize
  Superclass::Rasterize( mesh );

  // Make sure everything has gone smoothly
  if ( m_Abort )
    {
    // Something has gone wrong
    m_MinLogLikelihoodTimesPrior = itk::NumericTraits< double >::max();
    return;
    }
    
  // Collect the results of all the threads
  for ( std::vector< double >::const_iterator  it = m_ThreadSpecificMinLogLikelihoodTimesPriors.begin();
        it != m_ThreadSpecificMinLogLikelihoodTimesPriors.end(); ++it )
    {
    m_MinLogLikelihoodTimesPrior += *it;
    }  
  
  m_PositionGradient = AtlasPositionGradientContainerType::New();
  AtlasPositionGradientType  zeroEntry( 0.0f );
  for ( AtlasMesh::PointsContainer::ConstIterator pointIt = mesh->GetPoints()->Begin();
        pointIt != mesh->GetPoints()->End(); ++pointIt )
    {
    m_PositionGradient->InsertElement( pointIt.Index(), zeroEntry );
    }
    
  for ( std::vector< AtlasPositionGradientContainerType::Pointer >::const_iterator  
             it = m_ThreadSpecificPositionGradients.begin();
        it != m_ThreadSpecificPositionGradients.end(); ++it )
    {
    AtlasPositionGradientContainerType::Iterator  sourceIt = ( *it )->Begin(); 
    AtlasPositionGradientContainerType::Iterator  targetIt = m_PositionGradient->Begin(); 
    for ( ; targetIt != m_PositionGradient->End(); ++sourceIt, ++targetIt )
      {
      targetIt.Value() += sourceIt.Value();
      } 
      
    } // End loop over all threads  
    
  
  // Alter the raw gradient to impose sliding boundary conditions
  this->ImposeSlidingBoundaryConditions( mesh );
    
}    
  
    
  
 

//
//
//
void
AtlasMeshToIntensityImageGradientCalculatorCPU
::SetMeshToImageTransform( const TransformType* meshToImageTransform )
{
  
  TransformType::Pointer  inScopeHolder = 0;
  if ( !meshToImageTransform )
    {
    //meshToImageTransform = TransformType::New();
    inScopeHolder = TransformType::New();
    meshToImageTransform = inScopeHolder;
    }
  
  // The sliding boundary conditions are implemented by multiplying the raw gradient with an appropriate
  // 3x3 matrix, mapping it into the allowable subspace. Since a vertex can move or not in x, y, and z
  // directions, there are 8 distinct matrices to perform the mapping. These are precomputed here.
  //
  // The index [0...7] is computed as follows: canMoveX*4 + canMoveY*2 + canMoveX
  //
  // The appropriate mapping matrix for each condition is given by a least-squares fit, resulting in
  // 
  //     correctionMatrix = C * inv( C' * C ) * C'
  //
  // where C has appropriate columns of the 3x3 matrix part of the meshToImageTransform (we don't care about
  // the translation part): if a certain direction cannot move, that column is missing from C
  //
  // Special cases are index 0 (cannot move in any of the three directions), in which case correctionMatrix=0,
  // and index 7 (can move freely in all three directions), in which case correctionMatrix=1. It's stupid to
  // actually explicitly multiply these correction matrices with the raw gradient, but I'm computing them here
  // anyway to avoid future programming errors.
  //
  SlidingBoundaryCorrectionMatrixType  correctionMatrix;
  
  // canMoveX=0, canMoveY=0, canMoveZ=0
  correctionMatrix.Fill( 0.0 );
  m_SlidingBoundaryCorrectionMatrices[ 0 ] = correctionMatrix;
  
  // canMoveX=0, canMoveY=0, canMoveZ=1
  vnl_matrix_fixed< double, 3, 1 >  C1;
  C1( 0, 0 ) = meshToImageTransform->GetMatrix()( 0, 2 );
  C1( 1, 0 ) = meshToImageTransform->GetMatrix()( 1, 2 );
  C1( 2, 0 ) = meshToImageTransform->GetMatrix()( 2, 2 );
  m_SlidingBoundaryCorrectionMatrices[ 1 ] = ( C1 * vnl_inverse( C1.transpose() * C1 ) * C1.transpose() );
  
  // canMoveX=0, canMoveY=1, canMoveZ=0
  vnl_matrix_fixed< double, 3, 1 >  C2;
  C2( 0, 0 ) = meshToImageTransform->GetMatrix()( 0, 1 );
  C2( 1, 0 ) = meshToImageTransform->GetMatrix()( 1, 1 );
  C2( 2, 0 ) = meshToImageTransform->GetMatrix()( 2, 1 );
  m_SlidingBoundaryCorrectionMatrices[ 2 ] = ( C2 * vnl_inverse( C2.transpose() * C2 ) * C2.transpose() );

  // canMoveX=0, canMoveY=1, canMoveZ=1
  vnl_matrix_fixed< double, 3, 2 >  C3;
  C3( 0, 0 ) = meshToImageTransform->GetMatrix()( 0, 1 );
  C3( 1, 0 ) = meshToImageTransform->GetMatrix()( 1, 1 );
  C3( 2, 0 ) = meshToImageTransform->GetMatrix()( 2, 1 );
  C3( 0, 1 ) = meshToImageTransform->GetMatrix()( 0, 2 );
  C3( 1, 1 ) = meshToImageTransform->GetMatrix()( 1, 2 );
  C3( 2, 1 ) = meshToImageTransform->GetMatrix()( 2, 2 );
  m_SlidingBoundaryCorrectionMatrices[ 3 ] = ( C3 * vnl_inverse( C3.transpose() * C3 ) * C3.transpose() );

  // canMoveX=1, canMoveY=0, canMoveZ=0
  vnl_matrix_fixed< double, 3, 1 >  C4;
  C4( 0, 0 ) = meshToImageTransform->GetMatrix()( 0, 0);
  C4( 1, 0 ) = meshToImageTransform->GetMatrix()( 1, 0 );
  C4( 2, 0 ) = meshToImageTransform->GetMatrix()( 2, 0 );
  m_SlidingBoundaryCorrectionMatrices[ 4 ] = ( C4 * vnl_inverse( C4.transpose() * C4 ) * C4.transpose() );

  // canMoveX=1, canMoveY=0, canMoveZ=1
  vnl_matrix_fixed< double, 3, 2 >  C5;
  C5( 0, 0 ) = meshToImageTransform->GetMatrix()( 0, 0 );
  C5( 1, 0 ) = meshToImageTransform->GetMatrix()( 1, 0 );
  C5( 2, 0 ) = meshToImageTransform->GetMatrix()( 2, 0 );
  C5( 0, 1 ) = meshToImageTransform->GetMatrix()( 0, 2 );
  C5( 1, 1 ) = meshToImageTransform->GetMatrix()( 1, 2 );
  C5( 2, 1 ) = meshToImageTransform->GetMatrix()( 2, 2 );
  m_SlidingBoundaryCorrectionMatrices[ 5 ] = ( C5 * vnl_inverse( C5.transpose() * C5 ) * C5.transpose() );

  // canMoveX=1, canMoveY=1, canMoveZ=0
  vnl_matrix_fixed< double, 3, 2 >  C6;
  C6( 0, 0 ) = meshToImageTransform->GetMatrix()( 0, 0 );
  C6( 1, 0 ) = meshToImageTransform->GetMatrix()( 1, 0 );
  C6( 2, 0 ) = meshToImageTransform->GetMatrix()( 2, 0 );
  C6( 0, 1 ) = meshToImageTransform->GetMatrix()( 0, 1 );
  C6( 1, 1 ) = meshToImageTransform->GetMatrix()( 1, 1 );
  C6( 2, 1 ) = meshToImageTransform->GetMatrix()( 2, 1 );
  m_SlidingBoundaryCorrectionMatrices[ 6 ] = ( C6 * vnl_inverse( C6.transpose() * C6 ) * C6.transpose() );

  // canMoveX=1, canMoveY=1, canMoveZ=1
  correctionMatrix.SetIdentity();
  m_SlidingBoundaryCorrectionMatrices[ 7 ] = correctionMatrix;
  
}

//
//
//
bool
AtlasMeshToIntensityImageGradientCalculatorCPU
::RasterizeTetrahedron( const AtlasMesh* mesh, 
                        AtlasMesh::CellIdentifier tetrahedronId,
                        int threadNumber )
{
  
  // If some other thread already encountered a folded tetrahedron, return
  // "false" immediately so that our own thread can also immediately return
  // in the superclass
  if ( m_Abort )
    {
    return false;
    }
  
  
  // We start with an empty gradient vector in each 
  AtlasPositionGradientType  gradientInVertex0( 0.0 );
  AtlasPositionGradientType  gradientInVertex1( 0.0 );
  AtlasPositionGradientType  gradientInVertex2( 0.0 );
  AtlasPositionGradientType  gradientInVertex3( 0.0 );
  double  priorPlusDataCost = 0.0;

  
  
  // ==========================================================
  //
  // Part I: Cache relevant elements of the vertices of this tetrahedron.
  // This includes alphas and position in each vertex; and information about
  // the reference position of the tetrahedron (needed for evaluating the
  // contribution of the prior
  //
  // ==========================================================
          
  //
  ReferenceTetrahedronInfo  info;
  mesh->GetCellData( tetrahedronId, &info );
 
  AtlasMesh::CellAutoPointer  cell;
  mesh->GetCell( tetrahedronId, cell );

  AtlasMesh::CellType::PointIdIterator  pit = cell->PointIdsBegin();
  const AtlasMesh::PointIdentifier  id0 = *pit;
  ++pit;
  const AtlasMesh::PointIdentifier  id1 = *pit;
  ++pit;
  const AtlasMesh::PointIdentifier  id2 = *pit;
  ++pit;
  const AtlasMesh::PointIdentifier  id3 = *pit;
  
  AtlasMesh::PointType p0;
  AtlasMesh::PointType p1;
  AtlasMesh::PointType p2;
  AtlasMesh::PointType p3;
  mesh->GetPoint( id0, &p0 );
  mesh->GetPoint( id1, &p1 );
  mesh->GetPoint( id2, &p2 );
  mesh->GetPoint( id3, &p3 );
  

  
  // ==========================================================
  //
  // Part II: Add contribution to cost and gradient of the deformation prior
  //
  // ==========================================================
  if ( !m_IgnoreDeformationPrior )
    {
    // Z is inv( [ p0 p1 p2 p3; 1 1 1 1 ] ) of the tetrahedron in reference position
    const double  referenceVolumeTimesK = info.m_ReferenceVolumeTimesK;

    const double  z11 = info.m_Z11;
    const double  z21 = info.m_Z21;
    const double  z31 = info.m_Z31;
    const double  z41 = info.m_Z41;

    const double  z12 = info.m_Z12;
    const double  z22 = info.m_Z22;
    const double  z32 = info.m_Z32;
    const double  z42 = info.m_Z42;

    const double  z13 = info.m_Z13;
    const double  z23 = info.m_Z23;
    const double  z33 = info.m_Z33;
    const double  z43 = info.m_Z43;


    //  Y = [ p0 p1 p2 p3; 1 1 1 1 ]
    const double y11 = p0[ 0 ];
    const double y21 = p0[ 1 ];
    const double y31 = p0[ 2 ];

    const double y12 = p1[ 0 ];
    const double y22 = p1[ 1 ];
    const double y32 = p1[ 2 ];

    const double y13 = p2[ 0 ];
    const double y23 = p2[ 1 ];
    const double y33 = p2[ 2 ];
    
    const double y14 = p3[ 0 ];
    const double y24 = p3[ 1 ];
    const double y34 = p3[ 2 ];

    
    //
    // Let's add Ashburner's prior cost for the tethrahedron deformation from its reference position
    //
    const double  m11 = z11*y11 + z21*y12 + z31*y13 + z41*y14;
    const double  m21 = z11*y21 + z21*y22 + z31*y23 + z41*y24;
    const double  m31 = z11*y31 + z21*y32 + z31*y33 + z41*y34;
    const double  m12 = z12*y11 + z22*y12 + z32*y13 + z42*y14;
    const double  m22 = z12*y21 + z22*y22 + z32*y23 + z42*y24;
    const double  m32 = z12*y31 + z22*y32 + z32*y33 + z42*y34;
    const double  m13 = z13*y11 + z23*y12 + z33*y13 + z43*y14;
    const double  m23 = z13*y21 + z23*y22 + z33*y23 + z43*y24;
    const double  m33 = z13*y31 + z23*y32 + z33*y33 + z43*y34;

    const double  detJ = m11 * ( m22*m33 - m32*m23 ) - m12 * ( m21*m33 - m31*m23 ) + m13 * ( m21*m32 - m31*m22 );
    if ( detJ <= 0 )
      {
      m_Abort = true;
      return false;
      }
    


    // Let's define K as inv( J ) * det( J )
    const double  k11 = ( m22*m33 - m23*m32 );
    const double  k12 = -( m12*m33 - m32*m13 );
    const double  k13 = ( m12*m23 - m22*m13 );
    const double  k21 = -( m21*m33 - m31*m23 );
    const double  k22 = ( m11*m33 - m13*m31 );
    const double  k23 = -( m11*m23 - m21*m13 );
    const double  k31 = ( m21*m32 - m31*m22 );
    const double  k32 = -( m11*m32 - m31*m12 );
    const double  k33 = ( m11*m22 - m12*m21 );

    // Trace of J' * J is actually the sum of the squares of the singular values of J: s1^2 + s2^2 + s3^2
    const double  sumOfSquaresOfSingularValuesOfJ = m11*m11 + m12*m12 + m13*m13 + m21*m21 + m22*m22 + m23*m23 + m31*m31 + m32*m32 + m33*m33;

    // Trace of ( inv(J) )' * inv( J ) is actually the sum of the squares of the reciprocals of the singular values
    // of J: 1/s1^2 + 1/s2^2 + 1/s3^2
    const double  traceOfKTransposeTimesK = ( k11*k11 + k12*k12 + k13*k13 + k21*k21 + k22*k22 + k23*k23 + k31*k31 + k32*k32 + k33*k33 );
    const double  sumOfSquaresOfReciprocalsOfSingularValuesOfJ = traceOfKTransposeTimesK / ( detJ * detJ );

    const double  priorCost = referenceVolumeTimesK * ( 1 + detJ ) *
                            ( sumOfSquaresOfSingularValuesOfJ + sumOfSquaresOfReciprocalsOfSingularValuesOfJ - 6 );
    priorPlusDataCost += priorCost;



    //
    // OK, now add contribution to derivatives of Ashburner's prior cost in each of the tetrahedron's vertices
    //
    const double  ddetJdm11 = m22*m33 - m32*m23;
    const double  ddetJdm21 = m13*m32 - m12*m33;
    const double  ddetJdm31 = m12*m23 - m13*m22;
    const double  ddetJdm12 = m31*m23 - m21*m33;
    const double  ddetJdm22 = m11*m33 - m13*m31;
    const double  ddetJdm32 = m13*m21 - m11*m23;
    const double  ddetJdm13 = m21*m32 - m31*m22;
    const double  ddetJdm23 = m12*m31 - m11*m32;
    const double  ddetJdm33 = m11*m22 - m12*m21;


    const double  tmp1 = referenceVolumeTimesK * ( sumOfSquaresOfSingularValuesOfJ + sumOfSquaresOfReciprocalsOfSingularValuesOfJ - 6 );
    const double  tmp2 = 2 * referenceVolumeTimesK * ( 1 + detJ );
    const double  tmp3 = pow( detJ,  3 );


    const double  dcostdm11 = tmp1 * ddetJdm11 + tmp2 * ( m11 + ( ( k22*m33 - k23*m23 - k32*m32 + k33*m22 ) * detJ - traceOfKTransposeTimesK * ddetJdm11 ) / tmp3 );
    const double  dcostdm21 = tmp1 * ddetJdm21 + tmp2 * ( m21 + ( ( -k21*m33 + k23*m13 + k31*m32 - k33*m12 ) * detJ - traceOfKTransposeTimesK * ddetJdm21 ) / tmp3 );
    const double  dcostdm31 = tmp1 * ddetJdm31 + tmp2 * ( m31 + ( ( k21*m23 - k22*m13 - k31*m22 + k32*m12 ) * detJ - traceOfKTransposeTimesK * ddetJdm31 ) / tmp3 );

    const double  dcostdm12 = tmp1 * ddetJdm12 + tmp2 * ( m12 + ( ( -k12*m33 + k13*m23 + k32*m31 - k33*m21 ) * detJ - traceOfKTransposeTimesK * ddetJdm12 ) / tmp3 );
    const double  dcostdm22 = tmp1 * ddetJdm22 + tmp2 * ( m22 + ( ( k11*m33 - k13*m13 - k31*m31 + k33*m11 ) * detJ - traceOfKTransposeTimesK * ddetJdm22 ) / tmp3 );
    const double  dcostdm32 = tmp1 * ddetJdm32 + tmp2 * ( m32 + ( ( -k11*m23 + k12*m13 + k31*m21 - k32*m11 ) * detJ - traceOfKTransposeTimesK * ddetJdm32 ) / tmp3 );

    const double  dcostdm13 = tmp1 * ddetJdm13 + tmp2 * ( m13 + ( ( k12*m32 - k13*m22 - k22*m31 + k23*m21 ) * detJ - traceOfKTransposeTimesK * ddetJdm13 ) / tmp3 );
    const double  dcostdm23 = tmp1 * ddetJdm23 + tmp2 * ( m23 + ( ( -k11*m32 + k13*m12 + k21*m31 -k23*m11 ) * detJ - traceOfKTransposeTimesK * ddetJdm23 ) / tmp3 );
    const double  dcostdm33 = tmp1 * ddetJdm33 + tmp2 * ( m33 + ( ( k11*m22 - k12*m12 - k21*m21 + k22*m11 ) * detJ - traceOfKTransposeTimesK * ddetJdm33 ) / tmp3 );


    const double  dcostdy11 = dcostdm11*z11 + dcostdm12*z12 + dcostdm13*z13;
    const double  dcostdy21 = dcostdm21*z11 + dcostdm22*z12 + dcostdm23*z13;
    const double  dcostdy31 = dcostdm31*z11 + dcostdm32*z12 + dcostdm33*z13;

    const double  dcostdy12 = dcostdm11*z21 + dcostdm12*z22 + dcostdm13*z23;
    const double  dcostdy22 = dcostdm21*z21 + dcostdm22*z22 + dcostdm23*z23;
    const double  dcostdy32 = dcostdm31*z21 + dcostdm32*z22 + dcostdm33*z23;

    const double  dcostdy13 = dcostdm11*z31 + dcostdm12*z32 + dcostdm13*z33;
    const double  dcostdy23 = dcostdm21*z31 + dcostdm22*z32 + dcostdm23*z33;
    const double  dcostdy33 = dcostdm31*z31 + dcostdm32*z32 + dcostdm33*z33;

    const double  dcostdy14 = dcostdm11*z41 + dcostdm12*z42 + dcostdm13*z43;
    const double  dcostdy24 = dcostdm21*z41 + dcostdm22*z42 + dcostdm23*z43;
    const double  dcostdy34 = dcostdm31*z41 + dcostdm32*z42 + dcostdm33*z43;


    // Add the stuff to the existing gradients in each of the tetrahedron's four vertices
    gradientInVertex0[ 0 ] += dcostdy11;
    gradientInVertex0[ 1 ] += dcostdy21;
    gradientInVertex0[ 2 ] += dcostdy31;
    
    gradientInVertex1[ 0 ] += dcostdy12;
    gradientInVertex1[ 1 ] += dcostdy22;
    gradientInVertex1[ 2 ] += dcostdy32;

    gradientInVertex2[ 0 ] += dcostdy13;
    gradientInVertex2[ 1 ] += dcostdy23;
    gradientInVertex2[ 2 ] += dcostdy33;

    gradientInVertex3[ 0 ] += dcostdy14;
    gradientInVertex3[ 1 ] += dcostdy24;
    gradientInVertex3[ 2 ] += dcostdy34;
     
    }

  
  
  // ==========================================================
  //
  // Part III: Add contribution to cost and gradient of the data term
  //
  // ==========================================================
  if( !m_OnlyDeformationPrior )
    {
    // Loop over all voxels within the tetrahedron and do The Right Thing  
    const AtlasAlphasType&  alphasInVertex0 = mesh->GetPointData()->ElementAt( id0 ).m_Alphas;
    const AtlasAlphasType&  alphasInVertex1 = mesh->GetPointData()->ElementAt( id1 ).m_Alphas;
    const AtlasAlphasType&  alphasInVertex2 = mesh->GetPointData()->ElementAt( id2 ).m_Alphas;
    const AtlasAlphasType&  alphasInVertex3 = mesh->GetPointData()->ElementAt( id3 ).m_Alphas;
    const int  numberOfClasses = alphasInVertex0.Size();
    TetrahedronInteriorConstIterator< LikelihoodFilterType::OutputPixelType >  it( m_LikelihoodFilter->GetOutput(), p0, p1, p2, p3 );
    for ( unsigned int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
      it.AddExtraLoading( alphasInVertex0[ classNumber ], 
                          alphasInVertex1[ classNumber ], 
                          alphasInVertex2[ classNumber ], 
                          alphasInVertex3[ classNumber ] );
      }  
    for ( ; !it.IsAtEnd(); ++it )
      {
      // Skip voxels for which nothing is known
      if ( it.Value().Size() == 0 )
        {
        //std::cout << "Skipping: " << it.Value().Size() << std::endl;
        continue;
        }
        
      //
      double likelihood = 0.0;
      double  xGradientBasis = 0.0;
      double  yGradientBasis = 0.0;
      double  zGradientBasis = 0.0;
      for ( unsigned int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
        {
        // Get the Gaussian likelihood of this class at the intensity of this pixel
        const double gauss = it.Value()[ classNumber ];
          
        // Add contribution of the likelihood
        likelihood += gauss * it.GetExtraLoadingInterpolatedValue( classNumber );
        
        //
        xGradientBasis += gauss * it.GetExtraLoadingNextRowAddition( classNumber );
        yGradientBasis += gauss * it.GetExtraLoadingNextColumnAddition( classNumber );
        zGradientBasis += gauss * it.GetExtraLoadingNextSliceAddition( classNumber );
        } // End loop over all classes
        
        
      //  Add contribution to log-likelihood
      likelihood = likelihood + 1e-15; //dont want to divide by zero
      priorPlusDataCost -= log( likelihood );


      //
      xGradientBasis /= likelihood;
      yGradientBasis /= likelihood;
      zGradientBasis /= likelihood;

      // Add contribution to gradient in vertex 0
      gradientInVertex0[ 0 ] += xGradientBasis * it.GetPi0();
      gradientInVertex0[ 1 ] += yGradientBasis * it.GetPi0();
      gradientInVertex0[ 2 ] += zGradientBasis * it.GetPi0();

      // Add contribution to gradient in vertex 1
      gradientInVertex1[ 0 ] += xGradientBasis * it.GetPi1();
      gradientInVertex1[ 1 ] += yGradientBasis * it.GetPi1();
      gradientInVertex1[ 2 ] += zGradientBasis * it.GetPi1();
      
      // Add contribution to gradient in vertex 2
      gradientInVertex2[ 0 ] += xGradientBasis * it.GetPi2();
      gradientInVertex2[ 1 ] += yGradientBasis * it.GetPi2();
      gradientInVertex2[ 2 ] += zGradientBasis * it.GetPi2();
      
      // Add contribution to gradient in vertex 3
      gradientInVertex3[ 0 ] += xGradientBasis * it.GetPi3();
      gradientInVertex3[ 1 ] += yGradientBasis * it.GetPi3();
      gradientInVertex3[ 2 ] += zGradientBasis * it.GetPi3();
      
      
      } // End loop over all voxels within the tetrahedron

    } // End test if we need to include the data term
  
  
  
  
  // ==========================================================
  //
  // Part IV: Add contributions to the cost and gradient of this tetrahedron
  // 
  // ==========================================================
  m_ThreadSpecificMinLogLikelihoodTimesPriors[ threadNumber ] += priorPlusDataCost;

  m_ThreadSpecificPositionGradients[ threadNumber ]->ElementAt( id0 ) += gradientInVertex0;
  m_ThreadSpecificPositionGradients[ threadNumber ]->ElementAt( id1 ) += gradientInVertex1;
  m_ThreadSpecificPositionGradients[ threadNumber ]->ElementAt( id2 ) += gradientInVertex2;
  m_ThreadSpecificPositionGradients[ threadNumber ]->ElementAt( id3 ) += gradientInVertex3;

  return true;
}





//
//
//
void
AtlasMeshToIntensityImageGradientCalculatorCPU
::ImposeSlidingBoundaryConditions( const AtlasMesh* mesh )
{

  AtlasPositionGradientContainerType::Iterator  gradientIt = m_PositionGradient->Begin();
  AtlasMesh::PointDataContainer::ConstIterator  paramIt = mesh->GetPointData()->Begin();
  for ( ; gradientIt != m_PositionGradient->End(); ++gradientIt, ++paramIt )
    {
    // Construct the index [0...7] of the 8 possible boundary conditions
    int  index = 0;
    if ( paramIt.Value().m_CanMoveX )
      {
      index += 4;
      }
    if ( paramIt.Value().m_CanMoveY )
      {
      index += 2;
      }
    if ( paramIt.Value().m_CanMoveZ )
      {
      index += 1;
      }
      
    // Multiply by the correct 3x3 matrix imposing the appropriate boundary conditions.
    // For index 7, this is always the identity matrix, so we'll just skip the explicit
    // multiplication there (is vast majority of vertices anyway, so better be fast there)
    if ( index < 7 )
      {
      gradientIt.Value() = m_SlidingBoundaryCorrectionMatrices[ index ] * gradientIt.Value();
      }
      
    } // End loop over all points in the mesh

  
}



} // end namespace kvl
