#include "kvlAtlasMeshPositionCostAndGradientCalculator.h"

#include <itkMath.h>
#include "vnl/vnl_matrix_fixed.h"
#include "vnl/vnl_inverse.h"


namespace kvl
{

//
//
//
AtlasMeshPositionCostAndGradientCalculator
::AtlasMeshPositionCostAndGradientCalculator()
{

  m_MinLogLikelihoodTimesPrior = 0;
  m_IgnoreDeformationPrior = false;
  m_OnlyDeformationPrior = false;
  m_PositionGradient = 0;
  m_Abort = false;
  m_BoundaryCondition = SLIDING;

  this->SetMeshToImageTransform( TransformType::New() );
  
}


//
//
//
AtlasMeshPositionCostAndGradientCalculator
::~AtlasMeshPositionCostAndGradientCalculator()
{
}




//
//
//
void
AtlasMeshPositionCostAndGradientCalculator
::Rasterize( const AtlasMesh* mesh )
{

  // Initialize from a clean slate
  m_Abort = false;
  m_PositionGradient = 0;
  m_MinLogLikelihoodTimesPrior = 0;

#if KVL_ENABLE_TIME_PROBE  
  itk::TimeProbe clock;
  clock.Start();
#endif
  
  
  bool allocateNewMemory = true;
  if ( m_ThreadSpecificPositionGradients.size() == this->GetNumberOfThreads() )
    {
    if ( m_ThreadSpecificPositionGradients[0]->Size() == mesh->GetPoints()->Size() )
      {
      allocateNewMemory = false;  
      }
    }
  
  
  
  
  if ( allocateNewMemory )
    {
    m_ThreadSpecificPositionGradients.clear();
    m_ThreadSpecificMinLogLikelihoodTimesPriors.clear();
      
    // For each thread, create an empty gradient and cost so that
    // different threads never interfere with one another
    for ( int threadNumber = 0; threadNumber < this->GetNumberOfThreads(); threadNumber++ )
      {
      // Initialize cost to zero for this thread
      m_ThreadSpecificMinLogLikelihoodTimesPriors.push_back( 0.0 );  
        
      // Create a container to hold the position gradient of this thread, and initialize to zero
      AtlasPositionGradientThreadAccumContainerType::Pointer  positionGradient = AtlasPositionGradientThreadAccumContainerType::New();
      AtlasPositionGradientThreadAccumType  zeroEntry( 0.0f );
      for ( AtlasMesh::PointsContainer::ConstIterator pointIt = mesh->GetPoints()->Begin();
            pointIt != mesh->GetPoints()->End(); ++pointIt )
        {
        positionGradient->InsertElement( pointIt.Index(), zeroEntry );
        }
      m_ThreadSpecificPositionGradients.push_back( positionGradient );
      
      } // End loop over threads
    }
  else
    {
    // Simply zero out existing memory  
    for ( int threadNumber = 0; threadNumber < this->GetNumberOfThreads(); threadNumber++ )
      {
      m_ThreadSpecificMinLogLikelihoodTimesPriors[ threadNumber ] = 0.0;
        
      // 
      AtlasPositionGradientThreadAccumContainerType::Pointer  positionGradient
                                       = m_ThreadSpecificPositionGradients[ threadNumber ];
      AtlasPositionGradientThreadAccumType  zeroEntry( 0.0f );
      for ( auto gradIt = positionGradient->Begin();
            gradIt != positionGradient->End(); ++gradIt )
        {
        gradIt.Value() = zeroEntry;
        }
      
      } // End loop over threads
      
    }
      
    
#if KVL_ENABLE_TIME_PROBE      
  clock.Stop();
  std::cout << "Time taken by initialization: " << clock.GetMean() << std::endl;
#endif 
    
  // Now rasterize
#if KVL_ENABLE_TIME_PROBE  
  clock.Reset();
  clock.Start();
  m_ThreadSpecificDataTermRasterizationTimers = std::vector< itk::TimeProbe >( this->GetNumberOfThreads() );
  m_ThreadSpecificPriorTermRasterizationTimers = std::vector< itk::TimeProbe >( this->GetNumberOfThreads() );
  m_ThreadSpecificOtherRasterizationTimers = std::vector< itk::TimeProbe >( this->GetNumberOfThreads() );
#endif  
  Superclass::Rasterize( mesh );
#if KVL_ENABLE_TIME_PROBE  
  clock.Stop();
  std::cout << "Time taken by actual rasterization: " << clock.GetMean() << std::endl;
  double  dataTermRasterizationTime = 0.0;
  double  priorTermRasterizationTime = 0.0;
  double  otherRasterizationTime = 0.0;
  for ( int threadNumber = 0; threadNumber < this->GetNumberOfThreads(); threadNumber++ )
    {
    dataTermRasterizationTime += m_ThreadSpecificDataTermRasterizationTimers[ threadNumber ].GetTotal();
    priorTermRasterizationTime += m_ThreadSpecificPriorTermRasterizationTimers[ threadNumber ].GetTotal();
    otherRasterizationTime += m_ThreadSpecificOtherRasterizationTimers[ threadNumber ].GetTotal();
    }
  std::cout << "     dataTermRasterizationTime: " <<  dataTermRasterizationTime << std::endl; 
  std::cout << "     priorTermRasterizationTime: " <<  priorTermRasterizationTime << std::endl; 
  std::cout << "     otherRasterizationTime: " <<  otherRasterizationTime << std::endl; 
  
  clock.Reset();
  clock.Start();
#endif  
  
  
  if ( true )
    {
    // Initialize gradient-to-return to zero
    m_PositionGradient = AtlasPositionGradientContainerType::New();
    AtlasPositionGradientType  zeroEntry( 0.0f );
    for ( AtlasMesh::PointsContainer::ConstIterator pointIt = mesh->GetPoints()->Begin();
          pointIt != mesh->GetPoints()->End(); ++pointIt )
      {
      m_PositionGradient->InsertElement( pointIt.Index(), zeroEntry );
      }
    }
  else
    {
    // Initialize gradient-to-return to zero
    AtlasPositionGradientType  zeroEntry( 0.0f );
    for ( auto gradIt = m_PositionGradient->Begin();
          gradIt != m_PositionGradient->End(); ++gradIt )
      {
      gradIt.Value() = zeroEntry;
      }
    
    }

  // Make sure everything has gone smoothly
  if ( m_Abort )
    {
    // Something has gone wrong
    m_MinLogLikelihoodTimesPrior = itk::NumericTraits< double >::max();
    return;
    }
    
  // Collect MinLogLikelihoodTimesPrior across all threads
  ThreadAccumDataType totalThreadMinLogLikelihoodTimesPrior = 0;
  for ( int threadNumber = 0; threadNumber < this->GetNumberOfThreads(); threadNumber++ )
    {
    const double typedValue = double(m_ThreadSpecificMinLogLikelihoodTimesPriors[ threadNumber ]);
    if ( std::isnan( typedValue ) || std::isinf( typedValue ) )
      {
      // Something has gone wrong
      m_MinLogLikelihoodTimesPrior = itk::NumericTraits< double >::max();
      return;
      }
      
    totalThreadMinLogLikelihoodTimesPrior += m_ThreadSpecificMinLogLikelihoodTimesPriors[ threadNumber ];
    }

  // Copy accumulator value to final MinLogLikelihoodTimesPrior
  m_MinLogLikelihoodTimesPrior = totalThreadMinLogLikelihoodTimesPrior;

  // Accumulate PositionGradient across all threads
  for ( int threadNumber = 1; threadNumber < this->GetNumberOfThreads(); threadNumber++ )
    {
    AtlasPositionGradientThreadAccumContainerType::ConstIterator threadIt = m_ThreadSpecificPositionGradients[ threadNumber ]->Begin();
    AtlasPositionGradientThreadAccumContainerType::Iterator firstThreadIt = m_ThreadSpecificPositionGradients[ 0 ]->Begin();
    for ( ; firstThreadIt != m_ThreadSpecificPositionGradients[0]->End(); ++threadIt, ++firstThreadIt )
      {
      firstThreadIt.Value() += threadIt.Value();
      } 
    }

  // Copy accumulated values to final PositionGradients
  AtlasPositionGradientThreadAccumContainerType::Iterator firstThreadGradientIt = m_ThreadSpecificPositionGradients[ 0 ]->Begin(); 
  AtlasPositionGradientContainerType::Iterator finalGradientIt = m_PositionGradient->Begin(); 
  for ( ; finalGradientIt != m_PositionGradient->End(); ++firstThreadGradientIt, ++finalGradientIt )
    {
    finalGradientIt.Value() = firstThreadGradientIt.Value();
    }

#if KVL_ENABLE_TIME_PROBE  
  clock.Stop();
  std::cout << "Time taken by result collection: " << clock.GetMean() << std::endl;
  clock.Reset();
  clock.Start();
#endif
  
  //
  this->PostProcessCostAndGradient( mesh );
  
  // Take care of the desired boundary conditions
  this->ImposeBoundaryCondition( mesh );
  
#if KVL_ENABLE_TIME_PROBE  
  clock.Stop();
  std::cout << "Time taken by boundary condition imposition: " << clock.GetMean() << std::endl;
#endif   
    
}    
  
    
  
 

//
//
//
void
AtlasMeshPositionCostAndGradientCalculator
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
AtlasMeshPositionCostAndGradientCalculator
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

#if KVL_ENABLE_TIME_PROBE  
  m_ThreadSpecificOtherRasterizationTimers[ threadNumber ].Start();
#endif
  
    
#if 0  
  // We start with an empty gradient vector in each 
  AtlasPositionGradientType  gradientInVertex0( 0.0 );
  AtlasPositionGradientType  gradientInVertex1( 0.0 );
  AtlasPositionGradientType  gradientInVertex2( 0.0 );
  AtlasPositionGradientType  gradientInVertex3( 0.0 );
  double  priorPlusDataCost = 0.0;
  
  // Cache relevant things about the tetrahedron
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
  
#else
  
  // Cache relevant things about the tetrahedron
  //ReferenceTetrahedronInfo  info;
  //mesh->GetCellData( tetrahedronId, &info );
  // Implements internally mesh->GetCellData()->GetElementIfIndexExists(cellId, data);
  // More efficient is ReferenceTetrahedronInfo&  info = mesh->GetCellData()->ElementAt(ElementIdentifier) 
  const ReferenceTetrahedronInfo&  info = mesh->GetCellData()->ElementAt( tetrahedronId );
 
  //AtlasMesh::CellAutoPointer  cell;
  //mesh->GetCell( tetrahedronId, cell );
  //AtlasMesh::CellType::PointIdIterator  pit = cell->PointIdsBegin();
  // Implements internally: 
  //      CellType* cellptr = 0;
  //      this->GetCells()->GetElementIfIndexExists(cellId, &cellptr);
  //      cellPointer.TakeNoOwnership(cellptr);
  AtlasMesh::CellType::PointIdIterator  pit = mesh->GetCells()->ElementAt( tetrahedronId )->PointIdsBegin();
  const AtlasMesh::PointIdentifier  id0 = *pit;
  ++pit;
  const AtlasMesh::PointIdentifier  id1 = *pit;
  ++pit;
  const AtlasMesh::PointIdentifier  id2 = *pit;
  ++pit;
  const AtlasMesh::PointIdentifier  id3 = *pit;
  

  //AtlasMesh::PointType p0;
  //AtlasMesh::PointType p1;
  //AtlasMesh::PointType p2;
  //AtlasMesh::PointType p3;
  //mesh->GetPoint( id0, &p0 );
  //mesh->GetPoint( id1, &p1 );
  //mesh->GetPoint( id2, &p2 );
  //mesh->GetPoint( id3, &p3 );
  // Implements internally mesh->GetPoints()->GetElementIfIndexExists(ptId, point);
  // More efficient is AtlasMesh::PointType&  p0 = mesh->GetPoints()->ElementAt( id0 );
  const AtlasMesh::PointType&  p0 = mesh->GetPoints()->ElementAt( id0 );
  const AtlasMesh::PointType&  p1 = mesh->GetPoints()->ElementAt( id1 );
  const AtlasMesh::PointType&  p2 = mesh->GetPoints()->ElementAt( id2 );
  const AtlasMesh::PointType&  p3 = mesh->GetPoints()->ElementAt( id3 );
  
  ThreadAccumDataType&  priorPlusDataCost = m_ThreadSpecificMinLogLikelihoodTimesPriors[ threadNumber ];

  AtlasPositionGradientThreadAccumType&  gradientInVertex0 = m_ThreadSpecificPositionGradients[ threadNumber ]->ElementAt( id0 );
  AtlasPositionGradientThreadAccumType&  gradientInVertex1 = m_ThreadSpecificPositionGradients[ threadNumber ]->ElementAt( id1 );
  AtlasPositionGradientThreadAccumType&  gradientInVertex2 = m_ThreadSpecificPositionGradients[ threadNumber ]->ElementAt( id2 );
  AtlasPositionGradientThreadAccumType&  gradientInVertex3 = m_ThreadSpecificPositionGradients[ threadNumber ]->ElementAt( id3 );

  
#endif  
  
#if KVL_ENABLE_TIME_PROBE  
  m_ThreadSpecificOtherRasterizationTimers[ threadNumber ].Stop();
#endif  

  
  // Add contribution to cost and gradient of the deformation prior
  if ( !m_IgnoreDeformationPrior )
    {
#if KVL_ENABLE_TIME_PROBE  
    m_ThreadSpecificPriorTermRasterizationTimers[ threadNumber ].Start();
#endif    
    if ( !this->AddPriorContributionOfTetrahedron( p0, p1, p2, p3,
                                                   info, 
                                                   priorPlusDataCost,
                                                   gradientInVertex0, 
                                                   gradientInVertex1, 
                                                   gradientInVertex2, 
                                                   gradientInVertex3 ) )
      {
      m_Abort = true;
#if KVL_ENABLE_TIME_PROBE  
      m_ThreadSpecificPriorTermRasterizationTimers[ threadNumber ].Stop();
#endif      
      return false;
      }
#if KVL_ENABLE_TIME_PROBE  
    m_ThreadSpecificPriorTermRasterizationTimers[ threadNumber ].Stop();
#endif 
    } // End test if we need to include prior term

  
  
  // Add contribution to cost and gradient of the data term
  if( !m_OnlyDeformationPrior )
    {
#if KVL_ENABLE_TIME_PROBE  
    m_ThreadSpecificDataTermRasterizationTimers[ threadNumber ].Start();
#endif
    
    const AtlasAlphasType&  alphasInVertex0 = mesh->GetPointData()->ElementAt( id0 ).m_Alphas;
    const AtlasAlphasType&  alphasInVertex1 = mesh->GetPointData()->ElementAt( id1 ).m_Alphas;
    const AtlasAlphasType&  alphasInVertex2 = mesh->GetPointData()->ElementAt( id2 ).m_Alphas;
    const AtlasAlphasType&  alphasInVertex3 = mesh->GetPointData()->ElementAt( id3 ).m_Alphas;
  
    this->AddDataContributionOfTetrahedron( p0, p1, p2, p3,
                                            alphasInVertex0, 
                                            alphasInVertex1, 
                                            alphasInVertex2, 
                                            alphasInVertex3,
                                            priorPlusDataCost,
                                            gradientInVertex0, 
                                            gradientInVertex1, 
                                            gradientInVertex2, 
                                            gradientInVertex3 );
    
#if KVL_ENABLE_TIME_PROBE  
    m_ThreadSpecificDataTermRasterizationTimers[ threadNumber ].Stop();
#endif    
    } // End test if we need to include the data term
  

  return true;
}



//
//
//
bool
AtlasMeshPositionCostAndGradientCalculator
::AddPriorContributionOfTetrahedron( const AtlasMesh::PointType& p0,
                                     const AtlasMesh::PointType& p1,
                                     const AtlasMesh::PointType& p2,
                                     const AtlasMesh::PointType& p3,
                                     const ReferenceTetrahedronInfo& info,
                                     ThreadAccumDataType&  priorPlusDataCost,
                                     AtlasPositionGradientThreadAccumType&  gradientInVertex0,
                                     AtlasPositionGradientThreadAccumType&  gradientInVertex1,
                                     AtlasPositionGradientThreadAccumType&  gradientInVertex2,
                                     AtlasPositionGradientThreadAccumType&  gradientInVertex3 )
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
    
  return true;
}  



//
//
//
void
AtlasMeshPositionCostAndGradientCalculator
::ImposeBoundaryCondition( const AtlasMesh* mesh )
{

  switch( m_BoundaryCondition ) 
    {
    case SLIDING: 
      {
      //std::cout << "SLIDING" << std::endl;
      this->ImposeSlidingBoundaryConditions( mesh );      
      break;
      } 
    case AFFINE: 
      {
      //std::cout << "AFFINE" << std::endl;
      this->ImposeAffineBoundaryConditions( mesh );  
      break;
      } 
    case TRANSLATION: 
      {
      //std::cout << "TRANSLATION" << std::endl;
      this->ImposeTranslationBoundaryConditions( mesh );  
      break;
      } 
    default:
      {
      //std::cout << "NONE" << std::endl;
      break;
      }
    }
    
}

    
//
//
//
void
AtlasMeshPositionCostAndGradientCalculator
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



//
//
//
void
AtlasMeshPositionCostAndGradientCalculator
::ImposeAffineBoundaryConditions( const AtlasMesh* mesh )
{
  // Make sure the gradient obeys an affine transformation model. This is accomplished 
  // by using an implicit model according to which a new position 3D y = (y1, y2, y3)^T 
  // is obtained from the current 3D position x = (x1, x2, x3)^T using an affine transformation:
  //
  //   y1 = a11 * x1 + a12 * x2 + a13 * x3 + t1
  //   y2 = a21 * x2 + a22 * x2 + a23 * x3 + t2
  //   y3 = a31 * x2 + a32 * x2 + a33 * x3 + t3
  //
  // where t = (t1, t2, t3)^T is a translation and the "a" parameters are elements of
  // a 3x3 matrix A.
  //
  // Using the notation X to denote a ( numberOfPoints x 4 ) matrix with on the n-th row
  // the position ( x1 x2 x3 1 ) of the n-th point, and G the ( numberOfPoints x 3 ) gradient,
  // the projection of G(:,1) -- the first column of G, i.e., the "left-direction" direction 
  // parameterized by the first row of the system above -- into the affine subspace spanned 
  // by the vectors of X is given 
  // 
  //   Gaffine(:,1) = X * ( X^T * X )^{-1} * X^T * G(:,1)
  //
  // This can be seen by the fact that the projected gradient is given by the constrained form
  // Gaffine(:,1) = X * ( a11' a12' a13' t1' ) = X * p; requiring the difference between the
  // original and the projected  gradient vector to be perpendicular to basis vectors yields
  //
  //    ( X * p - G(:,1) )^T * X = 0
  //
  // which yields 
  //
  //    p = ( X^T * X )^{-1} * X^T * G(:,1)
  //
  // The same argument goes for the two other directions (rows in the system above)
  //
  // One last thing: we essentially get access to y, whereas we need x to construct X. I'm
  // assuming that the first-ever time this function is called, t=( 0 0 0)^T and A = identity,
  // so that x = y in that case; from this X is constructed and stored for future use.
  //
  // One truly last thing: instead of storing X and working with the rows of A and t as 
  // implicit parameters, I'm changing basis, essentially working with, for each direction,
  // 4x1 implicit parameter vector pnew = C * p, where C is the Cholesky decomposition of   
  // the matrix-to-be-inverted X^T * X. This is accomplished by using Z = X * C^{-1} as basis vectors 
  // instead (which span the same affine subspace), which is computed using essentially a 
  // hand-coded QR decomposition (Gram-Schmidt orthonormalization) -- this basis has the advantage 
  // that Z^T * Z = C^{-1}^T * X^T * X * C^{-1} = C^{-1}^T * C^T * C * C^{-1} = I,
  // so no need to invert matrices anymore
  
  // Let's first construct X and G
  const int  numberOfRows = mesh->GetPoints()->Size();
  vnl_matrix< double >  G( numberOfRows, 3 );
  AtlasPositionGradientContainerType::Iterator  gradientIt = m_PositionGradient->Begin();
  for ( int  rowNumber = 0; rowNumber < numberOfRows; ++rowNumber, ++gradientIt )
    {
    //
    G( rowNumber, 0 ) = gradientIt.Value()[ 0 ];
    G( rowNumber, 1 ) = gradientIt.Value()[ 1 ];
    G( rowNumber, 2 ) = gradientIt.Value()[ 2 ];
    }

  if ( m_AffineProjectionMatrix.empty() )
    {
    m_AffineProjectionMatrix = vnl_matrix< double >( numberOfRows, 4 );
    AtlasMesh::PointsContainer::ConstIterator pointIt = mesh->GetPoints()->Begin();
    for ( int  rowNumber = 0; rowNumber < numberOfRows; ++rowNumber, ++pointIt )
      {
      //
      m_AffineProjectionMatrix( rowNumber, 0 ) = 1.0;
      m_AffineProjectionMatrix( rowNumber, 1 ) = pointIt.Value()[ 0 ];  
      m_AffineProjectionMatrix( rowNumber, 2 ) = pointIt.Value()[ 1 ];  
      m_AffineProjectionMatrix( rowNumber, 3 ) = pointIt.Value()[ 2 ];  
      }
      
    // Gram-Schmidt Orthonormalization up-front, so that we don't have to worry about inverting 
    // the same 4x4 matrix over and over 
    m_AffineProjectionMatrix.set_column( 0, m_AffineProjectionMatrix.get_column( 0 ).normalize() );
    for ( int columnNumber = 1; columnNumber < 4; columnNumber++ )
      {
      vnl_vector< double >  column = m_AffineProjectionMatrix.get_column( columnNumber );
      for ( int previousColumnNumber = 0; previousColumnNumber < columnNumber; previousColumnNumber++ )
        {
        column -= inner_product( column, m_AffineProjectionMatrix.get_column( previousColumnNumber ) ) 
                  * m_AffineProjectionMatrix.get_column( previousColumnNumber );
        }
      m_AffineProjectionMatrix.set_column( columnNumber, column.normalize() );
      }
      
    //std::cout << m_AffineProjectionMatrix.transpose() * m_AffineProjectionMatrix << std::endl; 
    
    }  
  
  
  // Compute  conformedG = X * X^T * G
  const  vnl_matrix< double >  conformedG = m_AffineProjectionMatrix * ( m_AffineProjectionMatrix.transpose() * G );
  //const  vnl_matrix< double >  conformedG 
  //       = m_X * ( vnl_inverse( m_X.transpose() * m_X ) * ( m_X.transpose() * G ) );
  
  // Now copy the result back into our gradient representation
  gradientIt = m_PositionGradient->Begin();
  for ( int  rowNumber = 0; rowNumber < numberOfRows; ++rowNumber, ++gradientIt )
   {
   gradientIt.Value()[ 0 ] = conformedG( rowNumber, 0 );  
   gradientIt.Value()[ 1 ] = conformedG( rowNumber, 1 );  
   gradientIt.Value()[ 2 ] = conformedG( rowNumber, 2 );  
   }  
    
}




//
//
//
void
AtlasMeshPositionCostAndGradientCalculator
::ImposeTranslationBoundaryConditions( const AtlasMesh* mesh )
{
  // Same as for affine case, but just simpler (X is just a (normalized) column of ones)
  
  // Let's first construct X and G
  const int  numberOfRows = mesh->GetPoints()->Size();
  vnl_matrix< double >  X( numberOfRows, 1 );
  vnl_matrix< double >  G( numberOfRows, 3 );
  AtlasMesh::PointsContainer::ConstIterator pointIt = mesh->GetPoints()->Begin();
  AtlasPositionGradientContainerType::Iterator  gradientIt = m_PositionGradient->Begin();
  for ( int  rowNumber = 0; rowNumber < numberOfRows; ++rowNumber, ++pointIt, ++gradientIt )
   {
   //
   X( rowNumber, 0 ) = 1.0;  
  
   //
   G( rowNumber, 0 ) = gradientIt.Value()[ 0 ];
   G( rowNumber, 1 ) = gradientIt.Value()[ 1 ];
   G( rowNumber, 2 ) = gradientIt.Value()[ 2 ];
   }     
  X.normalize_columns(); // Effectively performs Gram-Schmidt Orthonormalization
  
  // Compute  conformedG = X * X^T * G
  const  vnl_matrix< double >  conformedG = X * ( X.transpose() * G );
  
  // Now copy the result back into our gradient representation
  gradientIt = m_PositionGradient->Begin();
  for ( int  rowNumber = 0; rowNumber < numberOfRows; ++rowNumber, ++gradientIt )
   {
   gradientIt.Value()[ 0 ] = conformedG( rowNumber, 0 );  
   gradientIt.Value()[ 1 ] = conformedG( rowNumber, 1 );  
   gradientIt.Value()[ 2 ] = conformedG( rowNumber, 2 );  
   }  
    
}


} // end namespace kvl
