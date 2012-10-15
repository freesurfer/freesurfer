/**
 * @file  kvlAtlasMeshHamiltonianPositionSampler.cxx
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:38 $
 *    $Revision: 1.3 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */
#include "kvlAtlasMeshHamiltonianPositionSampler.h"

#include "vnl/vnl_sample.h"
#include "kvlAtlasMeshToIntensityImageGradientCalculator.h"


namespace kvl
{


//
//
//
AtlasMeshHamiltonianPositionSampler
::AtlasMeshHamiltonianPositionSampler()
{
  m_MeshCollection = 0;
  m_Image = 0;
  m_NumberOfWarmupSweeps = 1000;
  m_NumberOfBetweenSamplesSweeps = 1000;
  m_NumberOfSamples = 10;
  m_InitialPositionNumber = -1;

  m_MaximalTrackingTime = 4.0f;
  m_TimeStep = 0.4f;
  m_Mass = 1.0f;

}



//
//
//
AtlasMeshHamiltonianPositionSampler
::~AtlasMeshHamiltonianPositionSampler()
{
}




//
//
//
void
AtlasMeshHamiltonianPositionSampler
::PrintSelf( std::ostream& os, itk::Indent indent ) const
{
}



//
//
//
void
AtlasMeshHamiltonianPositionSampler
::Reseed( int seed )
{
  vnl_sample_reseed( seed );
}



//
//
//
AtlasPositionGradientContainerType::Pointer
AtlasMeshHamiltonianPositionSampler
::GetGradient( const AtlasMesh::PointsContainer*  position, double& cost ) const
{

  // Create a full mesh
  AtlasMesh::Pointer  mesh = AtlasMesh::New();
  mesh->SetPoints( const_cast< AtlasMesh::PointsContainer* >( position ) );
  mesh->SetCells( m_MeshCollection->GetCells() );
  mesh->SetPointData( m_MeshCollection->GetPointParameters() );
  mesh->SetCellData( const_cast< AtlasMesh::CellDataContainer* >( m_MeshCollection->GetReferenceTetrahedronInfos() ) );


  // Convert means and variances to correct format
  int  numberOfClasses = m_Means.size();
  itk::Array< float >   means( numberOfClasses );
  itk::Array< float >   variances( numberOfClasses );
  for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
  {
    means[ classNumber ] = m_Means[ classNumber ];
    variances[ classNumber ] = m_Variances[ classNumber ];
  }


  // Calculate the gradient and the cost
  AtlasMeshToIntensityImageGradientCalculator::Pointer  gradientCalculator =
    AtlasMeshToIntensityImageGradientCalculator::New();
  gradientCalculator->SetLabelImage( m_Image );
  gradientCalculator->SetMeans( means );
  gradientCalculator->SetVariances( variances );
  //gradientCalculator->SetIgnoreDeformationPrior( useAffine );
  gradientCalculator->Rasterize( mesh );
  cost = gradientCalculator->GetMinLogLikelihoodTimesPrior();

  // Return the gradient
  return gradientCalculator->GetPositionGradient();

}



#if 1

//
//
//
AtlasMeshCollection::Pointer
AtlasMeshHamiltonianPositionSampler
::
GetSamples()
{

  // Sanity check on input
  if ( !m_MeshCollection )
  {
    itkExceptionMacro( << "No mesh collection set." );
  }
  if ( !m_Image )
  {
    itkExceptionMacro( << "Trying to sample from posterior but no image set." );
  }



  // Initial position is reference position if negative initial position number was given;
  // the corresponding position otherwise
  AtlasMesh::PointsContainer::ConstIterator  sourceIt
  = this->GetMeshCollection()->GetReferencePosition()->Begin();
  AtlasMesh::PointsContainer::ConstIterator  sourceItEnd
  = this->GetMeshCollection()->GetReferencePosition()->End();
  if ( m_InitialPositionNumber >= 0 )
  {
    sourceIt = this->GetMeshCollection()->GetPositions()[ m_InitialPositionNumber ]->Begin();
    sourceItEnd = this->GetMeshCollection()->GetPositions()[ m_InitialPositionNumber ]->End();
  }
  AtlasMesh::PointsContainer::Pointer  position = AtlasMesh::PointsContainer::New();
  for ( ; sourceIt != sourceItEnd ; ++sourceIt )
  {
    position->InsertElement( sourceIt.Index(), sourceIt.Value() );
  }


  // Allocate space for returning mesh collection containing samples
  AtlasMeshCollection::Pointer  samples = AtlasMeshCollection::New();
  samples->SetCells( m_MeshCollection->GetCells() );
  samples->SetPointParameters( m_MeshCollection->GetPointParameters() );
  samples->SetReferencePosition( m_MeshCollection->GetReferencePosition() );
  samples->SetK( m_MeshCollection->GetK() );
  std::vector< AtlasMesh::PointsContainer::Pointer >  samplePositions;
  samples->SetPositions( samplePositions );


  // Loop over MCMC sweeps
  double  cost;
  AtlasPositionGradientContainerType::Pointer  gradient = this->GetGradient( position, cost );
  if ( cost == itk::NumericTraits< double >::max() )
  {
    itkExceptionMacro( << "Initial position is invalid." );
  }
  for ( unsigned int sweepNumber = 0;
        sweepNumber < ( m_NumberOfWarmupSweeps + m_NumberOfBetweenSamplesSweeps * ( m_NumberOfSamples - 1 ) + 1 );
        sweepNumber++ )
  {

    // Sample from momentum
    AtlasPositionGradientContainerType::Pointer  momentum = AtlasPositionGradientContainerType::New();
    for ( AtlasMesh::PointDataContainer::ConstIterator  it = m_MeshCollection->GetPointParameters()->Begin();
          it != m_MeshCollection->GetPointParameters()->End(); ++it )
    {
      AtlasPositionGradientType  entry;
      if ( it.Value().m_CanMoveX )
      {
        entry[ 0 ] = sqrt( m_Mass ) * vnl_sample_normal( 0, 1 );
      }
      else
      {
        entry[ 0 ] = 0;
      }
      if ( it.Value().m_CanMoveY )
      {
        entry[ 1 ] = sqrt( m_Mass ) * vnl_sample_normal( 0, 1 );
      }
      else
      {
        entry[ 1 ] = 0;
      }
      if ( it.Value().m_CanMoveZ )
      {
        entry[ 2 ] = sqrt( m_Mass ) * vnl_sample_normal( 0, 1 );
      }
      else
      {
        entry[ 2 ] = 0;
      }
      momentum->InsertElement( it.Index(), entry );
    }


    // Calculate the hamiltonian
    double  hamiltonian = cost;
    for ( AtlasPositionGradientContainerType::ConstIterator  it = momentum->Begin();
          it != momentum->End(); ++it )
    {
      hamiltonian += ( pow( it.Value()[ 0 ], 2 ) + pow( it.Value()[ 1 ], 2 ) + pow( it.Value()[ 2 ], 2 ) ) / m_Mass / 2.0f;
    }


    // Initialize the trial position and trial gradient to the current values
    AtlasMesh::PointsContainer::Pointer  trialPosition = AtlasMesh::PointsContainer::New();
    for ( AtlasMesh::PointsContainer::ConstIterator  sourceIt = position->Begin();
          sourceIt != position->End();
          ++sourceIt )
    {
      trialPosition->InsertElement( sourceIt.Index(), sourceIt.Value() );
    }
    AtlasPositionGradientContainerType::Pointer  trialGradient =  AtlasPositionGradientContainerType::New();
    for ( AtlasPositionGradientContainerType::ConstIterator  sourceIt = gradient->Begin();
          sourceIt != gradient->End();
          ++sourceIt )
    {
      trialGradient->InsertElement( sourceIt.Index(), sourceIt.Value() );
    }


    // Decide on how long we're gonna track this movement
    const float  trackingTime = vnl_sample_uniform( 0, 1 ) * m_MaximalTrackingTime;
    double  trialCost = 0.0;
    std::cout << "trackingTime: "<< trackingTime << std::endl;
    std::cout << "m_TimeStep: "<< m_TimeStep << std::endl;
    for ( float time = 0; time < trackingTime; time += m_TimeStep )
    {
      // Make half-step in momentum, and make step in trial position
      AtlasMesh::PointsContainer::Iterator  trialPositionIt = trialPosition->Begin();
      AtlasPositionGradientContainerType::Iterator  momentumIt = momentum->Begin();
      AtlasPositionGradientContainerType::ConstIterator  trialGradientIt = trialGradient->Begin();
      for ( ; trialPositionIt != trialPosition->End(); ++trialPositionIt, ++momentumIt, ++trialGradientIt )
      {
        for ( int i = 0; i < 3; i++ )
        {
          if ( momentumIt.Value()[ i ] )  // Test if we can move in this direction
          {
            momentumIt.Value()[ i ] -= ( m_TimeStep / 2 ) * trialGradientIt.Value()[ i ];
            trialPositionIt.Value()[ i ] += m_TimeStep * momentumIt.Value()[ i ] / m_Mass;
          }
        }

      }


      // Calculate the gradient in the new trial position
      trialGradient = this->GetGradient( trialPosition, trialCost );

      // If we've manouvred ourselves into an impossible situation, don't go any further
      if ( trialCost == itk::NumericTraits< double >::max() )
      {
        break;
      }

      // Make half-step in momentum
      momentumIt = momentum->Begin();
      trialGradientIt = trialGradient->Begin();
      for ( ; momentumIt != momentum->End(); ++momentumIt, ++trialGradientIt )
      {
        for ( int i = 0; i < 3; i++ )
        {
          if ( momentumIt.Value()[ i ] )  // Test if we can move in this direction
          {
            momentumIt.Value()[ i ] -= ( m_TimeStep / 2 ) * trialGradientIt.Value()[ i ];
          }
        }
      }

      //std::cout << "." << std::flush;
#if 1
      {
        double  trialHamiltonian = trialCost;
        for ( AtlasPositionGradientContainerType::ConstIterator  it = momentum->Begin();
              it != momentum->End(); ++it )
        {
          trialHamiltonian += ( pow( it.Value()[ 0 ], 2 ) + pow( it.Value()[ 1 ], 2 ) + pow( it.Value()[ 2 ], 2 ) ) / m_Mass / 2.0f;
        }

        std::cout << "  trialHamiltonian - hamiltonian: " << trialHamiltonian - hamiltonian
                  << " (" << trialCost - cost << " + " << ( trialHamiltonian - trialCost ) - ( hamiltonian - cost ) << ")" << std::endl;

      }
#endif

    } // End movement tracking
    //std::cout << std::endl;

    // Calculate the value of the Hamoltanian of the current trial (should theoretically be the same as
    // the one before movement tracking)
    double  trialHamiltonian = trialCost;
    if ( trialHamiltonian != itk::NumericTraits< double >::max() )
    {
      for ( AtlasPositionGradientContainerType::ConstIterator  it = momentum->Begin();
            it != momentum->End(); ++it )
      {
        trialHamiltonian += ( pow( it.Value()[ 0 ], 2 ) + pow( it.Value()[ 1 ], 2 ) + pow( it.Value()[ 2 ], 2 ) ) / m_Mass / 2.0f;
      }

    }


    // Decide whether to accept
    const  double  differenceInHamiltonian = trialHamiltonian - hamiltonian;
    bool  accept = false;
    if ( differenceInHamiltonian < 0 )
    {
      accept = true;
    }
    else if ( vnl_sample_uniform( 0, 1 ) < exp( -differenceInHamiltonian ) )
    {
      accept = true;
    }
    std::cout << std::endl;
    std::cout << "-->     differenceInHamiltonian: " << differenceInHamiltonian << std::endl;
    std::cout << "                         accept: " << accept << std::endl;
    std::cout << std::endl;
    std::cout << "             difference in cost: " << trialCost - cost << std::endl;
    std::cout << "             difference in energy: " << ( trialHamiltonian - trialCost ) - ( hamiltonian - cost ) << std::endl;



    // Do the necessary book-keeping if move is accepted
    if (  accept )
    {
      // Just for the heck of it, calculate how much you've actually moved
      float  maximalDeformation = 0.0f;
      AtlasMesh::PointsContainer::ConstIterator  posIt = position->Begin();
      AtlasMesh::PointsContainer::ConstIterator  trialPosIt = trialPosition->Begin();
      for ( ; posIt != position->End(); ++posIt, ++trialPosIt )
      {
        const float  deformation = ( trialPosIt.Value() - posIt.Value() ).GetNorm();

        if ( deformation > maximalDeformation )
        {
          maximalDeformation = deformation;
        }
      }
      std::cout << "             maximalDeformation: " << maximalDeformation << std::endl;


      //
      gradient = trialGradient;
      position = trialPosition;
      cost = trialCost;
    }


    // Test if we should memorize this sample
    if ( ( sweepNumber >= m_NumberOfWarmupSweeps ) &&
         !( ( sweepNumber - m_NumberOfWarmupSweeps ) % m_NumberOfBetweenSamplesSweeps ) )
    {
      int  sampleNumber = static_cast< int >( ( sweepNumber - m_NumberOfWarmupSweeps ) /
                                              m_NumberOfBetweenSamplesSweeps );
      std::cout << "sampleNumber: " << sampleNumber << " (sweepNumber: " << sweepNumber << ")" << std::endl;

      // Copy the current position to the samples container
      AtlasMesh::PointsContainer::Pointer  target = AtlasMesh::PointsContainer::New();
      for ( AtlasMesh::PointsContainer::ConstIterator  sourceIt = position->Begin();
            sourceIt != position->End();
            ++sourceIt )
      {
        target->InsertElement( sourceIt.Index(), sourceIt.Value() );
      }
      samplePositions.push_back( target );
      samples->SetPositions( samplePositions );

#if 1
      {
        std::ostringstream  fileNameStream;
        fileNameStream << "intermediateSamples_" << sampleNumber+1 << "_samples.txt";
        samples->Write( fileNameStream.str().c_str() );
      }
#endif

    } // End save of current position


    std::cout << std::endl;
    std::cout << std::endl;

  } // End loop over all MCMC sweeps

  return samples;

}


#else

//
//
//
AtlasMeshCollection::Pointer
AtlasMeshHamiltonianPositionSampler
::
GetSamples()
{

  // Sanity check on input
  if ( !m_MeshCollection )
  {
    itkExceptionMacro( << "No mesh collection set." );
  }
  if ( !m_Image )
  {
    itkExceptionMacro( << "Trying to sample from posterior but no image set." );
  }



  // Initial position is reference position if negative initial position number was given;
  // the corresponding position otherwise
  AtlasMesh::PointsContainer::ConstIterator  sourceIt
  = this->GetMeshCollection()->GetReferencePosition()->Begin();
  AtlasMesh::PointsContainer::ConstIterator  sourceItEnd
  = this->GetMeshCollection()->GetReferencePosition()->End();
  if ( m_InitialPositionNumber >= 0 )
  {
    sourceIt = this->GetMeshCollection()->GetPositions()[ m_InitialPositionNumber ]->Begin();
    sourceItEnd = this->GetMeshCollection()->GetPositions()[ m_InitialPositionNumber ]->End();
  }
  AtlasMesh::PointsContainer::Pointer  position = AtlasMesh::PointsContainer::New();
  for ( ; sourceIt != sourceItEnd ; ++sourceIt )
  {
    position->InsertElement( sourceIt.Index(), sourceIt.Value() );
  }


  // Allocate space for returning mesh collection containing samples
  AtlasMeshCollection::Pointer  samples = AtlasMeshCollection::New();
  samples->SetCells( m_MeshCollection->GetCells() );
  samples->SetPointParameters( m_MeshCollection->GetPointParameters() );
  samples->SetReferencePosition( m_MeshCollection->GetReferencePosition() );
  samples->SetK( m_MeshCollection->GetK() );
  std::vector< AtlasMesh::PointsContainer::Pointer >  samplePositions;
  samples->SetPositions( samplePositions );


  // Loop over MCMC sweeps
  double  cost;
  AtlasPositionGradientContainerType::Pointer  gradient = this->GetGradient( position, cost );
  if ( cost == itk::NumericTraits< double >::max() )
  {
    itkExceptionMacro( << "Initial position is invalid." );
  }
  for ( unsigned int sweepNumber = 0;
        sweepNumber < ( m_NumberOfWarmupSweeps + m_NumberOfBetweenSamplesSweeps * ( m_NumberOfSamples - 1 ) + 1 );
        sweepNumber++ )
  {

    // Sample from momentum
    AtlasPositionGradientContainerType::Pointer  momentum = AtlasPositionGradientContainerType::New();
    for ( AtlasMesh::PointDataContainer::ConstIterator  it = m_MeshCollection->GetPointParameters()->Begin();
          it != m_MeshCollection->GetPointParameters()->End(); ++it )
    {
      AtlasPositionGradientType  entry;
#if 0
      {
        if ( it.Index() != 389 )
        {
          const_cast< bool& >( it.Value().m_CanMoveX ) = false;
          const_cast< bool& >( it.Value().m_CanMoveY ) = false;
          const_cast< bool& >( it.Value().m_CanMoveZ ) = false;
        }
      }
#endif
      if ( it.Value().m_CanMoveX )
      {
        entry[ 0 ] = sqrt( m_Mass ) * vnl_sample_normal( 0, 1 );
      }
      else
      {
        entry[ 0 ] = 0;
      }
      if ( it.Value().m_CanMoveY )
      {
        entry[ 1 ] = sqrt( m_Mass ) * vnl_sample_normal( 0, 1 );
      }
      else
      {
        entry[ 1 ] = 0;
      }
      if ( it.Value().m_CanMoveZ )
      {
        entry[ 2 ] = sqrt( m_Mass ) * vnl_sample_normal( 0, 1 );
      }
      else
      {
        entry[ 2 ] = 0;
      }
      momentum->InsertElement( it.Index(), entry );
    }


    // Calculate the hamiltonian
    double  hamiltonian = cost;
    for ( AtlasPositionGradientContainerType::ConstIterator  it = momentum->Begin();
          it != momentum->End(); ++it )
    {
      hamiltonian += ( pow( it.Value()[ 0 ], 2 ) + pow( it.Value()[ 1 ], 2 ) + pow( it.Value()[ 2 ], 2 ) ) / m_Mass / 2.0f;
    }


    // Initialize the trial position and trial gradient to the current values
    AtlasMesh::PointsContainer::Pointer  trialPosition = AtlasMesh::PointsContainer::New();
    for ( AtlasMesh::PointsContainer::ConstIterator  sourceIt = position->Begin();
          sourceIt != position->End();
          ++sourceIt )
    {
      trialPosition->InsertElement( sourceIt.Index(), sourceIt.Value() );
    }
    AtlasPositionGradientContainerType::Pointer  trialGradient =  AtlasPositionGradientContainerType::New();
    for ( AtlasPositionGradientContainerType::ConstIterator  sourceIt = gradient->Begin();
          sourceIt != gradient->End();
          ++sourceIt )
    {
      trialGradient->InsertElement( sourceIt.Index(), sourceIt.Value() );
    }


    // Decide on how long we're gonna track this movement
    const float  trackingTime = vnl_sample_uniform( 0, 1 ) * m_MaximalTrackingTime;
    double  trialCost = 0.0f;
    std::cout << "trackingTime: "<< trackingTime << std::endl;
    std::cout << "m_TimeStep: "<< m_TimeStep << std::endl;
    for ( float time = 0; time < trackingTime; time += m_TimeStep )
    {
      // Make a half-step in trial position
      AtlasMesh::PointsContainer::Iterator  trialPositionIt = trialPosition->Begin();
      AtlasPositionGradientContainerType::Iterator  momentumIt = momentum->Begin();
      for ( ; trialPositionIt != trialPosition->End(); ++trialPositionIt, ++momentumIt )
      {
#if 0
        if ( momentumIt.Index() == 389 )
        {
          std::cout << "trialPosition before: " << trialPositionIt.Value() << std::endl;
        }
#endif

        for ( int i = 0; i < 3; i++ )
        {
          if ( momentumIt.Value()[ i ] )  // Test if we can move in this direction
          {
            trialPositionIt.Value()[ i ] += ( m_TimeStep / 2.0f ) * momentumIt.Value()[ i ] / m_Mass;
          }
        }

#if 0
        if ( momentumIt.Index() == 389 )
        {
          std::cout << "trialPosition after: " << trialPositionIt.Value() << std::endl;
        }
#endif
      }

      // Calculate the gradient in the new trial position
      trialGradient = this->GetGradient( trialPosition, trialCost );

      // If we've manouvred ourselves into an impossible situation, don't go any further
      if ( trialCost == itk::NumericTraits< double >::max() )
      {
        break;
      }


      // Evaluate hypothetical hamiltonian of making only a half-step in trial momentum, multiplied
      // by a factor close to 1. Theoretically this should be equal to the original hamiltonian when
      // the factor is exactly 1
      float  bestAlpha = 1.0f;
      double  bestError = itk::NumericTraits< double >::max();
      const float  delta = 0.01;
      for ( float alpha = 0.9; alpha < 1.1; alpha += delta )
      {
        momentumIt = momentum->Begin();
        AtlasPositionGradientContainerType::ConstIterator  trialGradientIt = trialGradient->Begin();
        double  testHamiltonian = trialCost;
        for ( ; momentumIt != momentum->End(); ++momentumIt, ++trialGradientIt )
        {
          for ( int i = 0; i < 3; i++ )
          {
            if ( momentumIt.Value()[ i ] )  // Test if we can move in this direction
            {
              testHamiltonian += ( pow( ( momentumIt.Value()[ i ] - alpha * ( m_TimeStep / 2.0f ) * trialGradientIt.Value()[ i ] ), 2 ) / m_Mass / 2.0f );
            }
          }

        }

        //std::cout << "alpha: " << alpha << std::endl;
        //std::cout << "    testHamiltonian - hamiltonian: " << testHamiltonian - hamiltonian
        //          << " (" << trialCost - cost << " + " << ( testHamiltonian - trialCost ) - ( hamiltonian - cost ) << ")" << std::endl;

        const double  error = testHamiltonian - hamiltonian;
        if ( fabs( error ) < fabs( bestError ) )
        {
          bestAlpha = alpha;
          bestError = error;
        }
      } // End look for best alpha

      std::cout << "bestAlpha: " << bestAlpha << " (bestError: " << bestError << ")" << std::endl;


      // Make a step in trial momentum with the best alpha
      momentumIt = momentum->Begin();
      AtlasPositionGradientContainerType::ConstIterator  trialGradientIt = trialGradient->Begin();
      for ( ; momentumIt != momentum->End(); ++momentumIt, ++trialGradientIt )
      {
        for ( int i = 0; i < 3; i++ )
        {
          if ( momentumIt.Value()[ i ] )  // Test if we can move in this direction
          {
            momentumIt.Value()[ i ] -= bestAlpha * m_TimeStep * trialGradientIt.Value()[ i ];
          }
        }
      }


      // Make a half-step in trial position
      trialPositionIt = trialPosition->Begin();
      momentumIt = momentum->Begin();
      for ( ; trialPositionIt != trialPosition->End(); ++trialPositionIt, ++momentumIt )
      {
        for ( int i = 0; i < 3; i++ )
        {
          if ( momentumIt.Value()[ i ] )  // Test if we can move in this direction
          {
            trialPositionIt.Value()[ i ] += ( m_TimeStep / 2.0f ) * momentumIt.Value()[ i ] / m_Mass;
          }
        }

      }



      //std::cout << "." << std::flush;
    } // End movement tracking
    //std::cout << std::endl;

    // Calculate the value of the Hamoltanian of the current trial (should theoretically be the same as
    // the one before movement tracking)
    double  trialHamiltonian = trialCost;
    if ( trialHamiltonian != itk::NumericTraits< float >::max() )
    {
      for ( AtlasPositionGradientContainerType::ConstIterator  it = momentum->Begin();
            it != momentum->End(); ++it )
      {
        trialHamiltonian += ( pow( it.Value()[ 0 ], 2 ) + pow( it.Value()[ 1 ], 2 ) + pow( it.Value()[ 2 ], 2 ) ) / m_Mass / 2.0f;
      }

    }


    // Decide whether to accept
    const  double  differenceInHamiltonian = trialHamiltonian - hamiltonian;
    bool  accept = false;
    if ( differenceInHamiltonian < 0 )
    {
      accept = true;
    }
    else if ( vnl_sample_uniform( 0, 1 ) < exp( -differenceInHamiltonian ) )
    {
      accept = true;
    }
    std::cout << std::endl;
    std::cout << "-->     differenceInHamiltonian: " << differenceInHamiltonian << std::endl;
    std::cout << "                         accept: " << accept << std::endl;
    std::cout << std::endl;
    std::cout << "             difference in cost: " << trialCost - cost << std::endl;
    std::cout << "             difference in energy: " << ( trialHamiltonian - trialCost ) - ( hamiltonian - cost ) << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;



    // Do the necessary book-keeping if move is accepted
    if (  accept )
    {
      gradient = trialGradient;
      position = trialPosition;
      cost = trialCost;
    }


    // Test if we should memorize this sample
    if ( ( sweepNumber >= m_NumberOfWarmupSweeps ) &&
         !( ( sweepNumber - m_NumberOfWarmupSweeps ) % m_NumberOfBetweenSamplesSweeps ) )
    {
      int  sampleNumber = static_cast< int >( ( sweepNumber - m_NumberOfWarmupSweeps ) /
                                              m_NumberOfBetweenSamplesSweeps );
      std::cout << "sampleNumber: " << sampleNumber << " (sweepNumber: " << sweepNumber << ")" << std::endl;

      // Copy the current position to the samples container
      AtlasMesh::PointsContainer::Pointer  target = AtlasMesh::PointsContainer::New();
      for ( AtlasMesh::PointsContainer::ConstIterator  sourceIt = position->Begin();
            sourceIt != position->End();
            ++sourceIt )
      {
        target->InsertElement( sourceIt.Index(), sourceIt.Value() );
      }
      samplePositions.push_back( target );
      samples->SetPositions( samplePositions );

#if 1
      {
        std::ostringstream  fileNameStream;
        fileNameStream << "intermediateSamples_" << sampleNumber+1 << "_samples.txt";
        samples->Write( fileNameStream.str().c_str() );
      }
#endif

    } // End save of current position



  } // End loop over all MCMC sweeps

  return samples;

}

#endif

} // end namespace kvl
