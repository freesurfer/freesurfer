/**
 * @file  kvlAtlasMeshCollectionPositionCostCalculator.cxx
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
#include "kvlAtlasMeshCollectionPositionCostCalculator.h"



namespace kvl
{


//
//
//
AtlasMeshCollectionPositionCostCalculator
::AtlasMeshCollectionPositionCostCalculator()
{

}



//
//
//
AtlasMeshCollectionPositionCostCalculator
::~AtlasMeshCollectionPositionCostCalculator()
{

}



//
//
//
void
AtlasMeshCollectionPositionCostCalculator
::SetLabelImages( const std::vector< LabelImageType::ConstPointer >& labelImages )
{
  if ( labelImages.size() == 0 )
  {
    return;
  }

  m_LabelImages = labelImages;
  m_NumberOfLabelImages = m_LabelImages.size();

}







//
//
//
float
AtlasMeshCollectionPositionCostCalculator
::GetPositionCost( int& numberOfProblems )
{
  // Sanity check on input
  if ( ( !this->GetMeshCollection() ) || ( m_LabelImages.size() == 0 ) )
  {
    itkExceptionMacro( << "No mesh collection or label images set." );
  }

  // Let's check if we are expected to just return 0
  if ( this->GetMeshCollection()->GetK() >= 1000 )
  {
    return 0.0f;
  }

  // Initialize the number of problems to zero
  numberOfProblems = 0;

  // Loop over all images, adding their costs
  float  pooledPriorCurvatureCost = 0;
  float  pooledPosteriorCurvatureCost = 0;
  float  pooledPriorPositionCost = 0;
  float  pooledPosteriorPositionCost = 0;
  for ( unsigned int  labelImageNumber = 0; labelImageNumber < m_NumberOfLabelImages; labelImageNumber++ )
  {
    if ( this->GetDebug() )
    {
      std::cout << "Analyzing position cost for label image " << labelImageNumber << std::endl;
    }

    // Get curvature contributions of the data term
    PosteriorPositionCurvatureCalculator::Pointer  meshCalculator = PosteriorPositionCurvatureCalculator::New();
    meshCalculator->SetLabelImage( m_LabelImages[ labelImageNumber ] );
    meshCalculator->Rasterize( this->GetMeshCollection()->GetMesh( labelImageNumber ) );
    AtlasPositionCurvatureContainerType::ConstPointer  dataCurvatures = meshCalculator->GetPosteriorPositionCurvature();


    // Loop over all vertices, adding each vertex's contribution to the various cost components
    AtlasMesh::PointsContainer::Pointer  positions = this->GetMeshCollection()->GetPositions()[ labelImageNumber ];
    AtlasMesh::PointsContainer::ConstIterator  pointIt = positions->Begin();
    AtlasPositionCurvatureContainerType::ConstIterator  dataCurvIt = dataCurvatures->Begin();
    for ( ; pointIt != positions->End(); ++pointIt, ++dataCurvIt )
    {
      // Get the contribution of the prior MRF to the curvature in the current point
      const Curvature  priorCurvature = this->CalculateCurvature( labelImageNumber, pointIt.Index(),
                                        pointIt.Value()[ 0 ], pointIt.Value()[ 1 ] );

      // Check what components of the position are actually free parameters
      const bool  canMoveX = !( priorCurvature.m_Curvature_dxdx == itk::NumericTraits< float >::max() );
      const bool  canMoveY = !( priorCurvature.m_Curvature_dydy == itk::NumericTraits< float >::max() );
      if ( !canMoveX && !canMoveY )
      {
        // Nothing to do here
        continue;
      }

      // Get (xstar, ystar) and Cstar
      float xstar;
      float ystar;
      if ( !this->CalculateXstar( labelImageNumber, pointIt.Index(), xstar, ystar ) )
      {
        return itk::NumericTraits< float >::max();
      }
      const Curvature  CurvatureStar = this->CalculateCurvature( labelImageNumber, pointIt.Index(), xstar, ystar );
      const float  Cstarxx = CurvatureStar.m_Curvature_dxdx;
      const float  Cstarxy = CurvatureStar.m_Curvature_dxdy;
      const float  Cstaryy = CurvatureStar.m_Curvature_dydy;

      // Get the cost associated with the actual position of this vertex
      const float  posteriorPositionCost = this->CalculateCost( labelImageNumber, pointIt.Index(),
                                           pointIt.Value()[ 0 ], pointIt.Value()[ 1 ] );

      // Get the cost associated with the star position of this vertex
      const float  priorPositionCost = this->CalculateCost( labelImageNumber, pointIt.Index(),
                                       xstar, ystar );

      // Get the determinant of the Hessian for both the prior and the posterior
      float  detOfPosteriorCurvature = 0;
      float  detOfPriorCurvature = 0;
      if ( canMoveX && canMoveY )
      {
        // Normal case
        const float posteriorCurvature_dxdx = dataCurvIt.Value().m_Curvature_dxdx + priorCurvature.m_Curvature_dxdx;
        const float posteriorCurvature_dxdy = dataCurvIt.Value().m_Curvature_dxdy + priorCurvature.m_Curvature_dxdy;
        const float posteriorCurvature_dydy = dataCurvIt.Value().m_Curvature_dydy + priorCurvature.m_Curvature_dydy;
        detOfPosteriorCurvature = posteriorCurvature_dxdx * posteriorCurvature_dydy - pow( posteriorCurvature_dxdy, 2 );

        detOfPriorCurvature = Cstarxx * Cstaryy - Cstarxy * Cstarxy;
#if 0
        if ( detOfPosteriorCurvature < 0 )
        {
          // Try with the cross-term zeroed out
          detOfPosteriorCurvature = posteriorCurvature_dxdx * posteriorCurvature_dydy;
        }
#endif
      }
      else if ( canMoveX  )
      {
        // Only x can move
        detOfPosteriorCurvature = dataCurvIt.Value().m_Curvature_dxdx + priorCurvature.m_Curvature_dxdx;
        detOfPriorCurvature = Cstarxx;
      }
      else
      {
        // Only y can move
        detOfPosteriorCurvature = dataCurvIt.Value().m_Curvature_dydy + priorCurvature.m_Curvature_dydy;
        detOfPriorCurvature = Cstaryy;
      }


      // Calculate the actual cost components and this vertex, and it the total cost to the pooled cost
      pooledPriorPositionCost += priorPositionCost;
      pooledPosteriorPositionCost += posteriorPositionCost;
      pooledPriorCurvatureCost += 0.5 * log(  detOfPriorCurvature );
      if ( detOfPosteriorCurvature > 0 )
      {
        pooledPosteriorCurvatureCost += 0.5 * log( detOfPosteriorCurvature );
      }
      else
      {
        numberOfProblems++;
        //return itk::NumericTraits< float >::max();
      }


      // Print out some numbers
      if ( this->GetDebug() )
      {
        std::cout << "            point " << pointIt.Index() << ": " << std::endl;
        std::cout << "                 canMoveX: " << ( canMoveX ? "On" : "Off") << std::endl;
        std::cout << "                 canMoveY: " << ( canMoveY ? "On" : "Off") << std::endl;
        std::cout << "                 priorPositionCost: " << priorPositionCost << std::endl;
        std::cout << "                 posteriorPositionCost: " << posteriorPositionCost << std::endl;
        std::cout << "                 priorCurvatureCost: " <<  0.5 * log( detOfPriorCurvature ) << std::endl;
        std::cout << "                 posteriorCurvatureCost: " <<  0.5 * log( detOfPosteriorCurvature ) << std::endl;
      }


    } // End loop over all vertices

  } // End loop over all label images


  // Print out some collective info
  if ( this->GetDebug() )
  {
    std::cout << "=============================" << std::endl;
    std::cout << "Contribution to position cost: " << std::endl;
    std::cout << "       pooledPosteriorPositionCost: " << pooledPosteriorPositionCost << std::endl;
    std::cout << "           pooledPriorPositionCost: " << pooledPriorPositionCost << std::endl;
    std::cout << "      pooledPosteriorCurvatureCost: " << pooledPosteriorCurvatureCost << std::endl;
    std::cout << "          pooledPriorCurvatureCost: " << pooledPriorCurvatureCost << std::endl;
    std::cout << "  +/- -----------------------------------" << std::endl;
    std::cout << "                                    "
              << pooledPosteriorPositionCost + pooledPosteriorCurvatureCost -
              pooledPriorPositionCost - pooledPriorCurvatureCost << std::endl;
    std::cout << "=============================" << std::endl;
  }

  // Return result
  return ( pooledPosteriorPositionCost + pooledPosteriorCurvatureCost -
           pooledPriorPositionCost - pooledPriorCurvatureCost );

}




} // end namespace kvl
