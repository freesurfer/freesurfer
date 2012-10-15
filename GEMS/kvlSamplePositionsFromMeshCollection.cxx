/**
 * @file  kvlSamplePositionsFromMeshCollection.cxx
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
#include "kvlAtlasMeshHamiltonianPositionSampler.h"
#include "kvlEMSegmenter.h"
#include "itkImageFileReader.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkCastImageFilter.h"
#include <iostream>



int main( int argc, char** argv )
{
  // Sanity check
  if ( argc < 9 )
  {
    std::cerr << "Usage: " << argv[ 0 ]
              << " meshFileName imageFileName numberOfSamples numberOfWarmupSweeps numberOfBetweenSamplesSweeps mass"
              " timeStep maximalTrackingTime [initialPositionNumber seed]" << std::endl;
    return -1;
  }

  // Retrieve the input arguments
  std::ostringstream  inputParserStream;
  for ( int argumentNumber = 1; argumentNumber < argc; argumentNumber++ )
  {
    inputParserStream << argv[ argumentNumber ] << " ";
  }
  std::istringstream  inputStream( inputParserStream.str().c_str() );
  std::string  meshFileName;
  std::string  imageFileName;
  //unsigned int  domainSize[ 2 ];
  unsigned int  numberOfSamples;
  unsigned int  numberOfWarmupSweeps;
  unsigned int  numberOfBetweenSamplesSweeps;
  float  mass;
  float  timeStep;
  float  maximalTrackingTime;
  inputStream >>  meshFileName >> imageFileName >> numberOfSamples >> numberOfWarmupSweeps >> numberOfBetweenSamplesSweeps
              >> mass >> timeStep >> maximalTrackingTime;
  int  initialPositionNumber = -1;
  if ( argc > 9 )
  {
    inputStream >> initialPositionNumber;
    std::cout << "Explicitly set initialPositionNumber: " << initialPositionNumber << std::endl;
  }
  bool  reseed = false;
  int  seed = 0;
  if ( argc > 10 )
  {
    inputStream >> seed;
    std::cout << "Explicitly using seed " << seed << std::endl;
    reseed = true;
  }


  // Read the mesh collection
  kvl::AtlasMeshCollection::Pointer  collection = kvl::AtlasMeshCollection::New();
  collection->Read( meshFileName.c_str() );


  // Read the image
  typedef itk::Image< unsigned short, 3 >  ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  ReaderType::Pointer  reader = ReaderType::New();
  reader->SetFileName( imageFileName );
  reader->Update();
  ImageType::ConstPointer  image = reader->GetOutput();


  // Calculate uchar image
  std::cout << "TODO: we're using internally an image of unsigned char pixel type"
            " to calculate the position gradients" << std::endl;

  // Calculate minimum and maximum
  typedef itk::MinimumMaximumImageCalculator< ImageType >  RangeCalculatorType;
  RangeCalculatorType::Pointer  rangeCalculator = RangeCalculatorType::New();
  rangeCalculator->SetImage( image );
  rangeCalculator->Compute();
  std::cout << "minimum: " << rangeCalculator->GetMinimum() << std::endl;
  std::cout << "maximum: " << rangeCalculator->GetMaximum() << std::endl;

  // Scale and clip intensities to be between 0 and 255
  typedef kvl::AtlasMeshHamiltonianPositionSampler::ImageType  InternalImageType;
  typedef itk::IntensityWindowingImageFilter< ImageType, InternalImageType >   WindowerType;
  WindowerType::Pointer  windower = WindowerType::New();
  windower->SetInput( image );
  windower->SetWindowMinimum( rangeCalculator->GetMinimum() );
  windower->SetWindowMaximum( rangeCalculator->GetMaximum() );
  windower->SetOutputMinimum( 0 );
  windower->SetOutputMaximum( 255 );
  windower->Update();
  InternalImageType::ConstPointer  internalImage = windower->GetOutput();


  // Cast internal image back to original pixel type
  typedef itk::CastImageFilter< InternalImageType, ImageType >  BackCasterType;
  BackCasterType::Pointer  backCaster = BackCasterType::New();
  backCaster->SetInput( internalImage );
  backCaster->Update();
  image = backCaster->GetOutput();


  // Estimate means and variances
  kvl::EMSegmenter::Pointer  parameterEstimator = kvl::EMSegmenter::New();
  parameterEstimator->SetImage( image );
  if ( initialPositionNumber < 0 )
  {
    parameterEstimator->SetAtlasMesh( collection->GetReferenceMesh() );
  }
  else
  {
    parameterEstimator->SetAtlasMesh( collection->GetMesh( initialPositionNumber ) );
  }
  parameterEstimator->SetBiasFieldOrder( 0 );
  parameterEstimator->Segment();

  const std::vector< float >  means = parameterEstimator->GetMeans();
  const std::vector< float >  variances = parameterEstimator->GetVariances();

  std::cout << "means: [ ";
  for ( unsigned int i = 0; i < means.size(); i++ )
  {
    std::cout << means[ i ] << " ";
  }
  std::cout << "]" << std::endl;

  std::cout << "variances: [ ";
  for ( unsigned int i = 0; i < variances.size(); i++ )
  {
    std::cout << variances[ i ] << " ";
  }
  std::cout << "]" << std::endl;


  // Generate samples
  kvl::AtlasMeshHamiltonianPositionSampler::Pointer  sampler = kvl::AtlasMeshHamiltonianPositionSampler::New();
  sampler->SetMeshCollection( collection );
  sampler->SetImage( internalImage );
  sampler->SetMeans( means );
  sampler->SetVariances( variances );
  sampler->SetInitialPositionNumber( initialPositionNumber );
  sampler->SetNumberOfSamples( numberOfSamples );
  sampler->SetNumberOfWarmupSweeps( numberOfWarmupSweeps );
  sampler->SetNumberOfBetweenSamplesSweeps( numberOfBetweenSamplesSweeps );
  sampler->SetMass( mass );
  sampler->SetTimeStep( timeStep );
  sampler->SetMaximalTrackingTime( maximalTrackingTime );
  if ( reseed )
  {
    sampler->Reseed( seed );
  }
  kvl::AtlasMeshCollection::Pointer  samples = sampler->GetSamples();

  // Write the result out
  samples->Write( "positionSamples.txt" );

  return 0;
};


