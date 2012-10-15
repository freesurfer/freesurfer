/**
 * @file  kvlAtlasMeshCollectionPositionCostCalculator2.cxx
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
#include "kvlAtlasMeshCollectionPositionCostCalculator2.h"

#include "kvlAtlasMeshVertexProcessor2.h"
#include "vnl/vnl_det.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"


namespace kvl
{

bool AtlasMeshCollectionPositionCostCalculator::m_ReturnZero = false;


//
//
//
AtlasMeshCollectionPositionCostCalculator
::AtlasMeshCollectionPositionCostCalculator()
{

  m_MeshCollection = 0;

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
float
AtlasMeshCollectionPositionCostCalculator
::GetPositionCost( int& numberOfProblems )
{

  // Test if we're just supposed to return 0 quickly
  if ( m_ReturnZero )
  {
    numberOfProblems = 0;
    return 0.0f;
  }


  // Sanity check on input
  if ( !m_MeshCollection )
  {
    itkExceptionMacro( << "No mesh collection set." );
  }
  if ( m_MeshCollection->GetNumberOfMeshes() != m_LabelImages.size() )
  {
    itkExceptionMacro( << "Number of label images doesn't match number of meshes." );
  }

  // Let's check if we are expected to just return 0
  if ( m_MeshCollection->GetK() >= 1000 )
  {
    return 0.0f;
  }

  // Initiliaze the number of problems to zero
  numberOfProblems = 0;

  // Set up a vertex processor that will do all the work
  AtlasMeshVertexProcessor::Pointer  processor = AtlasMeshVertexProcessor::New();
  processor->SetMeshCollection( m_MeshCollection );
  processor->SetLabelImages( m_LabelImages );

  // Loop over all images, adding their costs
  float  pooledPriorCurvatureCost = 0;
  float  pooledPosteriorCurvatureCost = 0;
  float  pooledPriorPositionCost = 0;
  float  pooledPosteriorPositionCost = 0;
  for ( unsigned int  labelImageNumber = 0; labelImageNumber < m_LabelImages.size(); labelImageNumber++ )
  {
    if ( this->GetDebug() )
    {
      std::cout << "Analyzing position cost for label image " << labelImageNumber << std::endl;
    }

    // Loop over all vertices, adding each vertex's contribution to the various cost components
    AtlasMesh::PointsContainer::Pointer  positions = m_MeshCollection->GetPositions()[ labelImageNumber ];
    const int  numberOfNodes = positions->Size();
    int  numberOfNodesVisitedSoFar = 0;
    for ( AtlasMesh::PointsContainer::ConstIterator  pointIt = positions->Begin();
          pointIt != positions->End(); ++pointIt, ++numberOfNodesVisitedSoFar )
    {
      // Prior part
      processor->SetBeta( 0.0f );
      float  priorX;
      float  priorY;
      float  priorZ;
      if ( !processor->CalculateXstar( labelImageNumber, pointIt.Index(), priorX, priorY, priorZ ) )
      {
        //return itk::NumericTraits< float >::max();
      }
      const Curvature  priorCurvature = processor->CalculateCurvature( labelImageNumber, pointIt.Index(),
                                        priorX, priorY, priorZ );
      const float  priorCost = processor->CalculateCost( labelImageNumber, pointIt.Index(),
                               priorX, priorY, priorZ );


      // Posterior part
      processor->SetBeta( 1.0f );
#if 0
      const float  posteriorX = pointIt.Value()[ 0 ];
      const float  posteriorY = pointIt.Value()[ 1 ];
      const float  posteriorZ = pointIt.Value()[ 2 ];
#else
      float  posteriorX;
      float  posteriorY;
      float  posteriorZ;
      if ( !processor->CalculateXstar( labelImageNumber, pointIt.Index(), posteriorX, posteriorY, posteriorZ ) )
      {
        //return itk::NumericTraits< float >::max();
      }
#endif
      Curvature  posteriorCurvature = processor->CalculateCurvature( labelImageNumber, pointIt.Index(),
                                      posteriorX, posteriorY, posteriorZ );
      processor->SetBeta( 0.0f );
      const float  posteriorCost = processor->CalculateCost( labelImageNumber, pointIt.Index(),
                                   posteriorX, posteriorY, posteriorZ );

#if 1
      // Fiddle with the posterior curvature so that its determinant is positive. Following David MacKay's advice
      // we do this by seperating the data part from the deformation prior part, and ensuring that the data part
      // is positive semi-definite by setting any negative eigenvalues to zero. Unlike in MacKay's application,
      // where the prior part was guaranteed positive semi-definite because it was an isotropic Gaussian, our prior
      // part is *not* guaranteed to be, but typically it is.
      const Curvature  posteriorCurvatureContributionOfPrior =
        processor->CalculateCurvature( labelImageNumber, pointIt.Index(),
                                       posteriorX, posteriorY, posteriorZ );
      Curvature  posteriorCurvatureContributionOfData = posteriorCurvature;
      posteriorCurvatureContributionOfData -= posteriorCurvatureContributionOfPrior;


      // VNL stuff: first construct matrix, then calculate eigen decomposition, then change eigenvalues, and
      // finally recompose
      vnl_matrix< float >  posteriorCurvatureContributionOfDataMatrix( 3, 3 );
      posteriorCurvatureContributionOfDataMatrix( 0, 0 ) = posteriorCurvatureContributionOfData.m_Curvature_dxdx;
      posteriorCurvatureContributionOfDataMatrix( 0, 1 ) = posteriorCurvatureContributionOfData.m_Curvature_dxdy;
      posteriorCurvatureContributionOfDataMatrix( 0, 2 ) = posteriorCurvatureContributionOfData.m_Curvature_dxdz;
      posteriorCurvatureContributionOfDataMatrix( 1, 0 ) = posteriorCurvatureContributionOfData.m_Curvature_dxdy;
      posteriorCurvatureContributionOfDataMatrix( 1, 1 ) = posteriorCurvatureContributionOfData.m_Curvature_dydy;
      posteriorCurvatureContributionOfDataMatrix( 1, 2 ) = posteriorCurvatureContributionOfData.m_Curvature_dydz;
      posteriorCurvatureContributionOfDataMatrix( 2, 0 ) = posteriorCurvatureContributionOfData.m_Curvature_dxdz;
      posteriorCurvatureContributionOfDataMatrix( 2, 1 ) = posteriorCurvatureContributionOfData.m_Curvature_dydz;
      posteriorCurvatureContributionOfDataMatrix( 2, 2 ) = posteriorCurvatureContributionOfData.m_Curvature_dzdz;
      //std::cout << "original posteriorCurvatureContributionOfDataMatrix:\n"
      //          << posteriorCurvatureContributionOfDataMatrix << std::endl;

      vnl_symmetric_eigensystem< float >  eigenDecomposition( posteriorCurvatureContributionOfDataMatrix );
      for ( int i = 0 ; i < 3; i++ )
      {
        if ( eigenDecomposition.D( i ) < 0 )
        {
          std::cout << "Changing negative eigenvalue " << eigenDecomposition.D( i ) << " to zero" << std::endl;
          eigenDecomposition.D( i ) = 0;
        }
      }
      posteriorCurvatureContributionOfDataMatrix = eigenDecomposition.recompose();

      posteriorCurvatureContributionOfData.m_Curvature_dxdx = posteriorCurvatureContributionOfDataMatrix( 0, 0 );
      posteriorCurvatureContributionOfData.m_Curvature_dxdy = posteriorCurvatureContributionOfDataMatrix( 0, 1 );
      posteriorCurvatureContributionOfData.m_Curvature_dxdz = posteriorCurvatureContributionOfDataMatrix( 0, 2 );
      posteriorCurvatureContributionOfData.m_Curvature_dxdy = posteriorCurvatureContributionOfDataMatrix( 1, 0 );
      posteriorCurvatureContributionOfData.m_Curvature_dydy = posteriorCurvatureContributionOfDataMatrix( 1, 1 );
      posteriorCurvatureContributionOfData.m_Curvature_dydz = posteriorCurvatureContributionOfDataMatrix( 1, 2 );
      posteriorCurvatureContributionOfData.m_Curvature_dxdz = posteriorCurvatureContributionOfDataMatrix( 2, 0 );
      posteriorCurvatureContributionOfData.m_Curvature_dydz = posteriorCurvatureContributionOfDataMatrix( 2, 1 );
      posteriorCurvatureContributionOfData.m_Curvature_dzdz = posteriorCurvatureContributionOfDataMatrix( 2, 2 );
      //std::cout << "corrected posteriorCurvatureContributionOfDataMatrix:\n"
      //          << posteriorCurvatureContributionOfDataMatrix << std::endl;


      // Now recompose the posteriorCurvature from its parts
      posteriorCurvature = posteriorCurvatureContributionOfPrior;
      posteriorCurvature += posteriorCurvatureContributionOfData;

#endif

      // Get the determinant of the Hessian for both the prior and the posterior. We here make abuse of
      // the fact that points that are immobile automatically have their rows and columns filled with zeros,
      // except for the diagonals which have one. So we can just calculate the determinant and get it over
      // with
      vnl_matrix_fixed< float, 3, 3 >  priorCurvatureMatrix;
      priorCurvatureMatrix( 0, 0 ) = priorCurvature.m_Curvature_dxdx;
      priorCurvatureMatrix( 0, 1 ) = priorCurvature.m_Curvature_dxdy;
      priorCurvatureMatrix( 0, 2 ) = priorCurvature.m_Curvature_dxdz;
      priorCurvatureMatrix( 1, 0 ) = priorCurvature.m_Curvature_dxdy;
      priorCurvatureMatrix( 1, 1 ) = priorCurvature.m_Curvature_dydy;
      priorCurvatureMatrix( 1, 2 ) = priorCurvature.m_Curvature_dydz;
      priorCurvatureMatrix( 2, 0 ) = priorCurvature.m_Curvature_dxdz;
      priorCurvatureMatrix( 2, 1 ) = priorCurvature.m_Curvature_dydz;
      priorCurvatureMatrix( 2, 2 ) = priorCurvature.m_Curvature_dzdz;

      vnl_matrix_fixed< float, 3, 3 >  posteriorCurvatureMatrix;
      posteriorCurvatureMatrix( 0, 0 ) = posteriorCurvature.m_Curvature_dxdx;
      posteriorCurvatureMatrix( 0, 1 ) = posteriorCurvature.m_Curvature_dxdy;
      posteriorCurvatureMatrix( 0, 2 ) = posteriorCurvature.m_Curvature_dxdz;
      posteriorCurvatureMatrix( 1, 0 ) = posteriorCurvature.m_Curvature_dxdy;
      posteriorCurvatureMatrix( 1, 1 ) = posteriorCurvature.m_Curvature_dydy;
      posteriorCurvatureMatrix( 1, 2 ) = posteriorCurvature.m_Curvature_dydz;
      posteriorCurvatureMatrix( 2, 0 ) = posteriorCurvature.m_Curvature_dxdz;
      posteriorCurvatureMatrix( 2, 1 ) = posteriorCurvature.m_Curvature_dydz;
      posteriorCurvatureMatrix( 2, 2 ) = posteriorCurvature.m_Curvature_dzdz;

      const float  detOfPosteriorCurvature = vnl_det( posteriorCurvatureMatrix );
      const float  detOfPriorCurvature = vnl_det( priorCurvatureMatrix );


      // Calculate the actual cost components and this vertex, and it the total cost to the pooled cost
      pooledPriorPositionCost += priorCost;
      pooledPosteriorPositionCost += posteriorCost;
      if ( ( detOfPriorCurvature > 0 ) && ( detOfPosteriorCurvature > 0 ) )
      {
        pooledPriorCurvatureCost += 0.5 * log(  detOfPriorCurvature );
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
        std::cout << "          Progress: " << static_cast< float >( numberOfNodesVisitedSoFar ) /
                  static_cast< float >( numberOfNodes ) * 100 << "% ";
        std::cout << "            point " << pointIt.Index() << ": " << pointIt.Value() << std::endl;
        std::cout << "                 priorCost: " << priorCost << std::endl;
        std::cout << "                 posteriorCost: " << posteriorCost << std::endl;
        //std::cout << "                 priorCurvatureMatrix: \n" << priorCurvatureMatrix << std::endl;
        //std::cout << "                 posteriorCurvatureMatrix:\n" << posteriorCurvatureMatrix << std::endl;
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
