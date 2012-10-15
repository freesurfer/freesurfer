/**
 * @file  kvlAtlasMeshVertexProcessorTest.cxx
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:39 $
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
#include "kvlAtlasMeshCollection.h"
#include "kvlAtlasMeshVertexProcessor.h"
#include "itkImageFileReader.h"



int main( int argc, char* argv[] )
{
  if ( argc < 2 )
  {
    std::cerr << argv[0] << " meshFileName [beta imageFileName0 imageFileName1 ...]" << std::endl;
    return -1;
  }

  // Read the collection from file
  kvl::AtlasMeshCollection::Pointer  collection =  kvl::AtlasMeshCollection::New();
  collection->Read( argv[ 1 ] );

  // Set up the vertex processor
  kvl::AtlasMeshVertexProcessor::Pointer  processor = kvl::AtlasMeshVertexProcessor::New();
  processor->SetMeshCollection( collection );

  // Read beta and label images, and give them to the vertex processor
  if ( argc > 2 )
  {
    // Read beta
    std::istringstream  betaStream( argv[ 2 ] );
    float  beta;
    betaStream >> beta;
    processor->SetBeta( beta );

    // Read label images and give them to the vertex processor
    typedef itk::Image< unsigned char, 2 >  LabelImageType;
    std::vector< LabelImageType::ConstPointer >  labelImages;
    for ( int meshNumber = 0; meshNumber < collection->GetNumberOfMeshes(); meshNumber++ )
    {
      typedef itk::ImageFileReader< LabelImageType >  ReaderType;
      ReaderType::Pointer  reader = ReaderType::New();
      reader->SetFileName( argv[ 3 + meshNumber ] );
      reader->Update();
      labelImages.push_back( reader->GetOutput() );
    }
    processor->SetLabelImages( labelImages );
  }



  // Loop over all meshes
  for ( int meshNumber = 0; meshNumber < collection->GetNumberOfMeshes(); meshNumber++ )
  {
    std::cout << "=============================" << std::endl;
    std::cout << "meshNumber: " << meshNumber << std::endl;
    std::cout << "=============================" << std::endl;

    // Loop over all nodes
    for ( kvl::AtlasMesh::PointsContainer::ConstIterator  pointIt = collection->GetPositions()[ meshNumber ]->Begin();
          pointIt != collection->GetPositions()[ meshNumber ]->End(); ++pointIt )
    {
      const  kvl::AtlasMesh::PointIdentifier  pointId =  pointIt.Index();
      const  float  x = pointIt.Value()[ 0 ];
      const  float  y = pointIt.Value()[ 1 ];


      // Calculate the gradient and curvature analytically
      const kvl::AtlasPositionGradientType  analyticalGradient =
        processor->CalculateGradient( meshNumber, pointId, x, y );
      const kvl::Curvature  analyticalCurvature =
        processor->CalculateCurvature( meshNumber, pointId, x, y );

      // Calculate the gradient and curvature numerically
      const float  delta = 0.5 /* 0.005 */;
      const float  cost = processor->CalculateCost( meshNumber, pointId, x, y );
      const float  costPlusDeltaX = processor->CalculateCost( meshNumber, pointId, x + delta, y );
      const float  costMinusDeltaX = processor->CalculateCost( meshNumber, pointId, x - delta, y );
      const float  costPlusDeltaY = processor->CalculateCost( meshNumber, pointId, x, y + delta );
      const float  costMinusDeltaY = processor->CalculateCost( meshNumber, pointId, x, y - delta );
      const float  costPlusDeltaXPlusDeltaY = processor->CalculateCost( meshNumber, pointId, x + delta, y + delta );
      const float  costMinusDeltaXPlusDeltaY = processor->CalculateCost( meshNumber, pointId, x - delta, y + delta );
      const float  costPlusDeltaXMinusDeltaY = processor->CalculateCost( meshNumber, pointId, x + delta, y - delta );
      const float  costMinusDeltaXMinusDeltaY = processor->CalculateCost( meshNumber, pointId, x - delta, y - delta );

      kvl::AtlasPositionGradientType  numericalGradient;
      numericalGradient[ 0 ] = ( ( costPlusDeltaX - cost ) + ( cost - costMinusDeltaX ) ) / 2 / delta;
      numericalGradient[ 1 ] = ( ( costPlusDeltaY - cost ) + ( cost - costMinusDeltaY ) ) / 2 / delta;

      kvl::Curvature  numericalCurvature;
      numericalCurvature.m_Curvature_dxdx = ( ( costPlusDeltaX - cost ) - ( cost - costMinusDeltaX ) ) / pow( delta, 2 );
      const float  curvaturedxdyForward =  ( ( ( ( costPlusDeltaXPlusDeltaY - costPlusDeltaY ) +
                                             ( costPlusDeltaY - costMinusDeltaXPlusDeltaY ) ) / 2 / delta ) -
                                             ( ( ( costPlusDeltaX - cost ) +
                                                 ( cost - costMinusDeltaX ) ) / 2 / delta ) ) / delta;
      const float  curvaturedxdyBackward =  ( ( ( ( costPlusDeltaX - cost ) +
                                              ( cost - costMinusDeltaX ) ) / 2 / delta ) -
                                              ( ( ( costPlusDeltaXMinusDeltaY - costMinusDeltaY ) +
                                                  ( costMinusDeltaY - costMinusDeltaXMinusDeltaY ) ) / 2 / delta ) ) / delta;
      numericalCurvature.m_Curvature_dxdy = ( curvaturedxdyForward + curvaturedxdyBackward ) / 2;
      numericalCurvature.m_Curvature_dydy = ( ( costPlusDeltaY - cost ) - ( cost - costMinusDeltaY ) ) / pow( delta, 2 );


      // Print out results
      std::cout << "       Vertex with pointId: " << pointId << std::endl;
      std::cout << "                                x: " << x << std::endl;
      std::cout << "                                y: " << y << std::endl;
      std::cout << "                             cost: " << cost << std::endl;
      std::cout << std::endl;
      std::cout << "               analyticalGradient: " << analyticalGradient << std::endl;
      std::cout << "                numericalGradient: " << numericalGradient << std::endl;
      std::cout << std::endl;
      std::cout << "              analyticalCurvature: " << analyticalCurvature.m_Curvature_dxdx << ", "
                << analyticalCurvature.m_Curvature_dxdy << ", "
                << analyticalCurvature.m_Curvature_dydy << std::endl;
      std::cout << "               numericalCurvature: " << numericalCurvature.m_Curvature_dxdx << ", "
                << numericalCurvature.m_Curvature_dxdy << ", "
                << numericalCurvature.m_Curvature_dydy << std::endl;


    } // End loop over all nodes

  } // End loop over all meshes





  return 0;
};

