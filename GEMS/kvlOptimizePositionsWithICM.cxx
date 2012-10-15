/**
 * @file  kvlOptimizePositionsWithICM.cxx
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:40 $
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
#include "kvlAtlasMeshVertexProcessor2.h"
#include "itkImageFileReader.h"
#include "kvlAtlasMeshCollectionModelLikelihoodCalculator.h"
#include "kvlAtlasMeshCollectionPositionCostCalculator2.h"
#include "kvlAtlasParameterEstimator.h"


int main( int argc, char* argv[] )
{
  if ( argc < 4 )
  {
    std::cerr << argv[0] << " meshFileName K imageFileName0 imageFileName1 ..." << std::endl;
    return -1;
  }

  // Read the collection from file
  kvl::AtlasMeshCollection::Pointer  collection =  kvl::AtlasMeshCollection::New();
  collection->Read( argv[ 1 ] );

  // Retrieve K from input
  std::istringstream  KStream( argv[ 2 ] );
  float  K;
  KStream >> K;

  // Read label images
  typedef itk::Image< unsigned char, 3 >  LabelImageType;
  std::vector< LabelImageType::ConstPointer >  labelImages;
  for ( unsigned int meshNumber = 0; meshNumber < collection->GetNumberOfMeshes(); meshNumber++ )
  {
    typedef itk::ImageFileReader< LabelImageType >  ReaderType;
    ReaderType::Pointer  reader = ReaderType::New();
    reader->SetFileName( argv[ 3 + meshNumber ] );
    reader->Update();
    labelImages.push_back( reader->GetOutput() );
  }

#if 0
  float  currentK = 0.1;
#else
  float  currentK = K;
#endif
  while ( true )
  {
    // Set the K
    std::cout << "=============================" << std::endl;
    std::cout << "Analyzing for K: " << currentK << std::endl;
    std::cout << "=============================" << std::endl;
    collection->SetK( currentK );

    // Loop over all iterations
    for ( int iterationNumber = 0; iterationNumber < 30; iterationNumber++ )
    {
      std::cout << "  iterationNumber: " << iterationNumber << std::endl;

      // Set up the vertex processor
      kvl::AtlasMeshVertexProcessor::Pointer  processor = kvl::AtlasMeshVertexProcessor::New();
      processor->SetMeshCollection( collection );
      processor->SetLabelImages( labelImages );
      processor->SetBeta( 1.0 );

      // Loop over all meshes
      for ( unsigned int meshNumber = 0; meshNumber < collection->GetNumberOfMeshes(); meshNumber++ )
      {
        std::cout << "      meshNumber: " << meshNumber << std::endl;

        float  maximumMagnitudeOfChange = 0;

        // Loop over all nodes
        const int  numberOfNodes = collection->GetPositions()[ meshNumber ]->Size();
        int  numberOfNodesVisitedSoFar = 0;
        for ( kvl::AtlasMesh::PointsContainer::Iterator  pointIt = collection->GetPositions()[ meshNumber ]->Begin();
              pointIt != collection->GetPositions()[ meshNumber ]->End(); ++pointIt, ++numberOfNodesVisitedSoFar )
        {
          std::cout << "          Progress: " << static_cast< float >( numberOfNodesVisitedSoFar ) /
                    static_cast< float >( numberOfNodes ) * 100 << "% ";

          // Calculate xstar and ystar
          float  xstar;
          float  ystar;
          float  zstar;
          if ( processor->CalculateXstar( meshNumber, pointIt.Index(), xstar, ystar, zstar ) )
          {
            kvl::AtlasMesh::PointType  updatedPoint;
            updatedPoint[ 0 ] = xstar;
            updatedPoint[ 1 ] = ystar;
            updatedPoint[ 2 ] = zstar;
            const float  magnitudeOfChange = pointIt.Value().EuclideanDistanceTo( updatedPoint );
            if ( magnitudeOfChange > maximumMagnitudeOfChange )
            {
              maximumMagnitudeOfChange = magnitudeOfChange;
            }

            std::cout << "          Changing position of point " <<  pointIt.Index() << " from " << pointIt.Value()
                      << " to " << updatedPoint << " (magnitude: " << magnitudeOfChange << ")" << std::endl;

            pointIt.Value() = updatedPoint;
          }
          else
          {
            std::cout << "          Ooops.... couldn't update position of point " << pointIt.Index() << std::endl;
          }

        } // End loop over all nodes

        std::cout << "Maximal magnitude of change was " << maximumMagnitudeOfChange << std::endl;

      } // End loop over all meshes

#if 0
      {
        // Write out the current mesh collection
        std::ostringstream  fileNameStream;
        fileNameStream << "debug__optimizedWithICM_iteration" << iterationNumber << ".txt";
        collection->Write( fileNameStream.str().c_str() );
      }
#endif

      // Update the alphas
      std::cout << "   Updating alphas" << std::endl;
      collection->SetK( 10000 );
      //collection->FlattenAlphas();
      kvl::AtlasParameterEstimator::Pointer  estimator = kvl::AtlasParameterEstimator::New();
      estimator->SetLabelImages( labelImages );
      estimator->SetInitialMeshCollection( collection );
      estimator->Estimate();
      collection = const_cast< kvl::AtlasMeshCollection* >( estimator->GetCurrentMeshCollection() );
      collection->SetK( currentK );

      if ( false )
      {
        // Calculate the data cost and the alpha cost
        kvl::AtlasMeshCollectionModelLikelihoodCalculator::Pointer  dataAndAlphaCostCalculator =
          kvl::AtlasMeshCollectionModelLikelihoodCalculator::New();
        dataAndAlphaCostCalculator->SetMeshCollection( collection );
        dataAndAlphaCostCalculator->SetLabelImages( labelImages );
        float  dataCost;
        float  alphasCost;
        dataAndAlphaCostCalculator->GetDataCostAndAlphasCost( dataCost, alphasCost );

        // Calculate the position cost
        kvl::AtlasMeshCollectionPositionCostCalculator::Pointer  positionCostCalculator =
          kvl::AtlasMeshCollectionPositionCostCalculator::New();
        positionCostCalculator->SetMeshCollection( collection );
        positionCostCalculator->SetLabelImages( labelImages );
        positionCostCalculator->DebugOn();
        int  numberOfProblems;
        const float positionCost = positionCostCalculator->GetPositionCost( numberOfProblems );

        // Output total cost
        std::cout << "   Number of problems: "<< numberOfProblems << std::endl;
        std::cout << "   Total cost: " << std::endl;
        std::cout << "                    dataCost: " << dataCost << std::endl;
        std::cout << "                  alphasCost: " << alphasCost << std::endl;
        std::cout << "                positionCost: " << positionCost << std::endl;
        std::cout << "     + ---------------------  : " << std::endl;
        std::cout << "                    " << dataCost + alphasCost + positionCost << std::endl;
        std::cout << std::endl;
      }


      // Write out the current mesh collection
      std::ostringstream  fileNameStream;
      fileNameStream << "optimizedWithICM_" << currentK << "_iteration" << iterationNumber << ".txt";
      collection->Write( fileNameStream.str().c_str() );

    } // End loop over all iterations


    // Adjust K if we're not there yet. Otherwise, stop
    if ( currentK == K )
    {
      break;
    }
    currentK /= 1.4;
    if ( currentK <= K )
    {
      currentK = K;
    }

  } // End loop over K-adjusting steps

  return 0;
};






