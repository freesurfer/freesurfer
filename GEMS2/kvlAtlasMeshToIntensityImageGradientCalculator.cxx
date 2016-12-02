#include "kvlAtlasMeshToIntensityImageGradientCalculator.h"


namespace kvl
{

//
//
//
AtlasMeshToIntensityImageGradientCalculator
::AtlasMeshToIntensityImageGradientCalculator()
{
  
  this->SetMeshToImageTransform( TransformType::New() );
  
}


//
//
//
AtlasMeshToIntensityImageGradientCalculator
::~AtlasMeshToIntensityImageGradientCalculator()
{
}


//
//
//
double
AtlasMeshToIntensityImageGradientCalculator
::GetMinLogLikelihoodTimesPrior() const
{
  double  minLogLikelihoodTimesPrior = 0;
  for ( std::vector< FragmentProcessorType >::const_iterator  it = this->GetFragmentProcessors().begin();
        it != this->GetFragmentProcessors().end(); ++it )
    {
    if ( it->GetMinLogLikelihoodTimesPrior() == itk::NumericTraits< double >::max() )
      {
      return itk::NumericTraits< double >::max();  
      }

    minLogLikelihoodTimesPrior += it->GetMinLogLikelihoodTimesPrior();
    }
    
  return minLogLikelihoodTimesPrior;
}


  
  
//
//
//
AtlasPositionGradientContainerType::Pointer
AtlasMeshToIntensityImageGradientCalculator
::GetPositionGradient()
{
  
  //AtlasPositionGradientContainerType::ConstPointer  gradient = this->GetFragmentProcessor().GetPositionGradient();

  // Retrieve the "raw" gradient, i.e., without the sliding boundary conditions applied. This requires
  // summing the contributions over all threads.
  AtlasPositionGradientContainerType::Pointer  gradient = AtlasPositionGradientContainerType::New();
  AtlasPositionGradientType  zeroEntry( 0.0f );
  for ( AtlasPositionGradientContainerType::ConstIterator  it = this->GetFragmentProcessor().GetPositionGradient()->Begin();
        it != this->GetFragmentProcessor().GetPositionGradient()->End(); ++it )
    {
    gradient->InsertElement( it.Index(), zeroEntry );
    }
  for ( std::vector< FragmentProcessorType >::const_iterator it = this->GetFragmentProcessors().begin();
        it != this->GetFragmentProcessors().end(); ++it )
    {
    //// std::cout << "Adding contribution of thread to gradient" << std::endl;  
    AtlasPositionGradientContainerType::ConstIterator  sourceIt = it->GetPositionGradient()->Begin();
    AtlasPositionGradientContainerType::Iterator  targetIt = gradient->Begin();
    for ( ; targetIt != gradient->End(); ++targetIt, ++sourceIt )
      {
      targetIt.Value() += sourceIt.Value();
      }
    
    } // End loop over all threads
  //// std::cout << "Done looping over threads" << std::endl;  


  // Now loop over all vertices, look up how constrained they are by the sliding boundary conditions, and do the Right Thing
  //// std::cout << "Appling sliding boundary conditions" << std::endl;  
  AtlasPositionGradientContainerType::Iterator  gradientIt = gradient->Begin();
  AtlasMesh::PointDataContainer::ConstIterator  paramIt = this->GetFragmentProcessor().GetMesh()->GetPointData()->Begin();
  for ( ; gradientIt != gradient->End(); ++gradientIt, ++paramIt )
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
      
    //// std::cout << "index: " << index << std::endl;
        
    // Multiply by the correct 3x3 matrix imposing the appropriate boundary conditions.
    // For index 7, this is always the identity matrix, so we'll just skip the explicit
    // multiplication there (is vast majority of vertices anyway, so better be fast there)
    if ( index < 7 )
      {
      //// std::cout << "Applying boundary condition with index: "<< index << std::endl;
      //// std::cout << "    " << "entry before: " << entry << std::endl;
      gradientIt.Value() = m_SlidingBoundaryCorrectionMatrices[ index ] * gradientIt.Value();
      //// std::cout << "    " << "entry after: " << entry << std::endl;
      }
    }
  //// std::cout << "Done appling sliding boundary conditions" << std::endl;  


  // Return result
  return gradient;

}




//
//
//
void
AtlasMeshToIntensityImageGradientCalculator
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

  
  //
  //for ( int i=0; i<8; i++ )
  //  {
  //  // std::cout << "m_SlidingBoundaryCorrectionMatrices[ " << i << " ]: " << m_SlidingBoundaryCorrectionMatrices[ i ] << std::endl;
  //  }
  
}


} // end namespace kvl
