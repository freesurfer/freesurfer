#include <boost/test/unit_test.hpp>

#include "kvlAtlasMeshCollection.h"
#include "kvlAtlasMeshAlphaDrawer.h"
#include "kvlAtlasMeshVisitCounter.h"
#include "kvlAtlasMeshToIntensityImageCostAndGradientCalculator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageFileWriter.h"
#include "itkTimeProbe.h"

#include "testfileloader.hpp"

// -----------------------------------------

BOOST_FIXTURE_TEST_SUITE( AtlasMeshRasterizorBasic, TestFileLoader )

BOOST_AUTO_TEST_CASE( Interpolation )
{
  // Set up a timer
  itk::TimeProbe clock;

  // Rasterize the mesh, simply linearly interpolating a probabilistic atlas across the
  // volume of each tetrahedron
  kvl::AtlasMeshAlphaDrawer::Pointer  alphaDrawer = kvl::AtlasMeshAlphaDrawer::New();
  alphaDrawer->SetRegions( image->GetLargestPossibleRegion() );
  alphaDrawer->SetClassNumber( 1 );
  clock.Start();
  alphaDrawer->Rasterize( mesh );
  clock.Stop();
  BOOST_TEST_MESSAGE( "Time taken by alpha drawer: " << clock.GetMean() );
  

  // Write out
  typedef itk::ImageFileWriter< kvl::AtlasMeshAlphaDrawer::ImageType >  AlphaWriterType;
  AlphaWriterType::Pointer  alphaWriter = AlphaWriterType::New();
  alphaWriter->SetFileName( "testAlpha.nii" );
  alphaWriter->SetInput( alphaDrawer->GetImage() );
  alphaWriter->Write();
  BOOST_TEST_CHECKPOINT( "File written out" );

  
  // Compare against a reference implementation
  kvl::AtlasMeshAlphaDrawer::Pointer  referenceAlphaDrawer = kvl::AtlasMeshAlphaDrawer::New();
  referenceAlphaDrawer->SetRegions( image->GetLargestPossibleRegion() );
  referenceAlphaDrawer->SetClassNumber( 1 );
  clock.Reset();
  clock.Start();
  referenceAlphaDrawer->Rasterize( mesh );
  clock.Stop();
  BOOST_TEST_MESSAGE( "Time taken by reference alpha drawer: " << clock.GetMean() );
  itk::ImageRegionConstIteratorWithIndex< kvl::AtlasMeshAlphaDrawer::ImageType >  
             alphaIt( alphaDrawer->GetImage(), 
                      alphaDrawer->GetImage()->GetBufferedRegion() );
  itk::ImageRegionConstIteratorWithIndex< kvl::AtlasMeshAlphaDrawer::ImageType >  
             referenceAlphaIt( referenceAlphaDrawer->GetImage(), 
                               referenceAlphaDrawer->GetImage()->GetBufferedRegion() );
  double maximumAlphaError = 0.0;           
  for ( ; !alphaIt.IsAtEnd(); ++alphaIt, ++referenceAlphaIt ) {
    const double  error = std::abs( alphaIt.Value() - referenceAlphaIt.Value() );
    if ( error > maximumAlphaError ) {
      maximumAlphaError = error;
    }  
  }  
  BOOST_CHECK_EQUAL( maximumAlphaError, 0 );
}


BOOST_AUTO_TEST_CASE( VoxelCount )
{
  // Set up a timer
  itk::TimeProbe clock;

  kvl::AtlasMeshVisitCounter::Pointer  visitCounter = kvl::AtlasMeshVisitCounter::New();
  visitCounter->SetRegions( image->GetLargestPossibleRegion() );
  clock.Start();
  visitCounter->Rasterize( mesh );
  clock.Stop();
  BOOST_TEST_MESSAGE( "Time taken by visit counter: " << clock.GetMean() );

  itk::ImageRegionConstIteratorWithIndex< kvl::AtlasMeshVisitCounter::ImageType >  
             it( visitCounter->GetImage(), 
                 visitCounter->GetImage()->GetBufferedRegion() );
  bool success = true;
  for( ; !it.IsAtEnd(); ++it ) {
    BOOST_TEST_CONTEXT( "Voxel Index: " << it.GetIndex() ) {
      BOOST_CHECK_LE( it.Value(), 1 );
    }
  }

  // Write out
  typedef itk::ImageFileWriter< kvl::AtlasMeshVisitCounter::ImageType >  CountWriterType;
  CountWriterType::Pointer  countWriter = CountWriterType::New();
  countWriter->SetFileName( "testCount.nii" );
  countWriter->SetInput( visitCounter->GetImage() );
  countWriter->Write();
}

BOOST_AUTO_TEST_CASE( DeformationGradients )
{
  // Set up a timer
  itk::TimeProbe clock;
  
  std::vector< double >  means;
  means.push_back( 2.5429e3 );
  means.push_back( 3.2005e3 );
  means.push_back( 3.9733e3 );
  means.push_back( 4.6842e3 );
  means.push_back( 4.7763e3 );
  means.push_back( 4.2508e3 );
  means.push_back( 4.4670e3 );
  means.push_back( 4.5531e3 );
  means.push_back( 3.8819e3 );
  means.push_back( 4.0985e3 );
  means.push_back( 4.3246e3 );
  means.push_back( 4.6027e3 );
  means.push_back( 4.6872e3 );
  means.push_back( 4.5976e3 );
  means.push_back( 4.6632e3 );
  means.push_back( 4.7341e3 );
  means.push_back( 4.7243e3 );
 
  std::vector< double >  precisions;
  precisions.push_back( 0.000001230847303 );
  precisions.push_back( 0.000001179516824 );
  precisions.push_back( 0.000008872600132 );
  precisions.push_back( 0.000252957561955 );
  precisions.push_back( 0.001088434034911 );
  precisions.push_back( 0.000029275401400 );
  precisions.push_back( 0.000081278899442 );
  precisions.push_back( 0.000144651419610 );
  precisions.push_back( 0.000008551322644 );
  precisions.push_back( 0.000016721659636 );
  precisions.push_back( 0.000056766547482 );
  precisions.push_back( 0.000253226418815 );
  precisions.push_back( 0.000644518692120 );
  precisions.push_back( 0.000322204283476 );
  precisions.push_back( 0.000672630053763 );
  precisions.push_back( 0.000145138180789 );
  precisions.push_back( 0.000731213454142 );

  // Convert means and precisions to correct format   
  const int  numberOfClasses = means.size();
  std::vector< vnl_vector< double > > means_( numberOfClasses, vnl_vector< double >( 1, 0.0 ) );
  std::vector< vnl_matrix< double > >  precisions_( numberOfClasses, vnl_matrix< double >( 1, 1, 0.0 ) );
  for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ ) {
    means_[ classNumber ][0] = means[ classNumber ];
    precisions_[ classNumber ][0][0] = precisions[ classNumber ];
  }
  
  // Add in mixtureWeights and numberOfGaussiansPerClass
  std::vector<double> mixtureWeights( numberOfClasses );
  std::vector<int> numberOfGaussiansPerClass( numberOfClasses );
  for( int i=0; i<numberOfClasses; i++ ) {
    mixtureWeights.at(i) = numberOfGaussiansPerClass.at(i) = 1;
  }

  // Convert precisions to variances
  std::vector< vnl_matrix<double> > variances( precisions_ );
  for( int i=0; i<numberOfClasses; i++ ) {
    auto curr = variances.at(i);
    if( curr.size() != 1 ) {
      throw std::runtime_error("Must have 1x1 precisions matrix");
    }
    curr(0,0) = 1 / curr(0,0);
    variances.at(i) = curr;
  }

  //
  kvl::AtlasMeshToIntensityImageCostAndGradientCalculator::Pointer  
      gradientCalculator = kvl::AtlasMeshToIntensityImageCostAndGradientCalculator::New();
  gradientCalculator->SetImages( std::vector< ImageType::ConstPointer >( 1, image.GetPointer() ) );
  gradientCalculator->SetParameters( means_, variances, mixtureWeights, numberOfGaussiansPerClass );
  clock.Reset();
  clock.Start();
  gradientCalculator->Rasterize( mesh );
  clock.Stop();
  BOOST_TEST_MESSAGE( "Time taken for first call of gradient calculator: " << clock.GetMean() );

  // Let's do the timing also for subsequent iterations. This should be faster because an internal
  // ITK filter doesn't require updating until new images/means/precisions are set
  clock.Reset();
  for ( int testRunNumber = 0; testRunNumber < 10; testRunNumber++ ) {
    clock.Start();
    gradientCalculator->Rasterize( mesh );
    clock.Stop();
    BOOST_TEST_MESSAGE("Average time taken for subsequent iterations of gradient calculator: " << clock.GetMean() );
  }
      
      
  // Compare against a reference implementation
  kvl::AtlasMeshToIntensityImageCostAndGradientCalculator::Pointer  
      referenceGradientCalculator = kvl::AtlasMeshToIntensityImageCostAndGradientCalculator::New();
  referenceGradientCalculator->SetImages( std::vector< ImageType::ConstPointer >( 1, const_cast< ImageType* >( image.GetPointer() ) ) );
  referenceGradientCalculator->SetParameters( means_, variances, mixtureWeights, numberOfGaussiansPerClass );
  clock.Reset();
  clock.Start();
  referenceGradientCalculator->Rasterize( mesh );
  //referenceGradientCalculator->SetNumberOfThreads( 2 );
  clock.Stop();
  BOOST_TEST_MESSAGE("Time taken for first call of reference gradient calculator: " << clock.GetMean() );
  
  // Let's do the timing also for subsequent iterations. This should be faster because an internal
  // ITK filter doesn't require updating until new images/means/precisions are set
  clock.Reset();
  for ( int testRunNumber = 0; testRunNumber < 10; testRunNumber++ ) {
    clock.Start();
    referenceGradientCalculator->Rasterize( mesh );
    clock.Stop();
    BOOST_TEST_MESSAGE( "Average time taken for subsequent iterations of reference gradient calculator: " << clock.GetMean() );
  }
      
      
  // Compare results
  const double  cost = gradientCalculator->GetMinLogLikelihoodTimesPrior();
  const double  referenceCost = referenceGradientCalculator->GetMinLogLikelihoodTimesPrior();
  BOOST_CHECK_EQUAL( referenceCost, cost );
  
  kvl::AtlasPositionGradientContainerType::ConstPointer  gradient = gradientCalculator->GetPositionGradient();
  kvl::AtlasPositionGradientContainerType::ConstPointer  referenceGradient = referenceGradientCalculator->GetPositionGradient();
  kvl::AtlasPositionGradientContainerType::ConstIterator  gradIt = gradient->Begin();
  kvl::AtlasPositionGradientContainerType::ConstIterator  refGradIt = referenceGradient->Begin();
  double  maximumGradientError = 0.0;
  for ( ; gradIt != gradient->End(); ++gradIt, ++refGradIt ) {
    double  error = 0.0; 
    for ( int i = 0; i < 3; i++ ) {
      error += std::abs( gradIt.Value()[ i ] - refGradIt.Value()[ i ] );  
    }
    if ( error > maximumGradientError ) {
      maximumGradientError = error;
    }  
  }
  BOOST_CHECK_EQUAL( maximumGradientError, 0 );
}

BOOST_AUTO_TEST_SUITE_END();
