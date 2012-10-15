/**
 * @file  kvlReduceNumberOfIntensityLevels.cxx
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
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkThresholdLabelerImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkScalarImageToHistogramGenerator.h"




int main( int argc, char** argv )
{
  // Sanity check on the input
  if ( argc < 3 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " imageFileName numberOfIntensityLevels [ tagAlongImageFileName ... ]" << std::endl;
    return -1;
  }

  // Collect input parameters
  const std::string  imageFileName( argv[ 1 ] );
  std::istringstream  numberOfIntensityLevelsStream( argv[ 2 ] );
  int  numberOfIntensityLevels;
  numberOfIntensityLevelsStream >> numberOfIntensityLevels;

  //
  try
  {
    // Read the original image in
    typedef itk::Image< unsigned short, 3 >  ImageType;
    typedef itk::ImageFileReader< ImageType >  ReaderType;
    ReaderType::Pointer  reader = ReaderType::New();
    reader->SetFileName( imageFileName.c_str() );
    reader->Update();


    // Calculate min and max (this should probably be replaced by some percentile, robust method)
    typedef itk::MinimumMaximumImageCalculator< ImageType >  RangeCalculatorType;
    RangeCalculatorType::Pointer  rangeCalculator = RangeCalculatorType::New();
    rangeCalculator->SetImage( reader->GetOutput() );
    rangeCalculator->Compute();
    std::cout << "Minimum intensity: " << rangeCalculator->GetMinimum() << std::endl;
    std::cout << "Maximum intensity: " << rangeCalculator->GetMaximum() << std::endl;


    // Calculate the histogram of the image
    const int  numberOfBins = 256;
    //const float  binWidth = static_cast< float >( rangeCalculator->GetMaximum() - rangeCalculator->GetMinimum() + 1 ) /
    //                        static_cast< float >( numberOfBins );
    typedef itk::Statistics::ScalarImageToHistogramGenerator< ImageType >   HistogramGeneratorType;
    HistogramGeneratorType::Pointer histogramGenerator = HistogramGeneratorType::New();
    histogramGenerator->SetInput(  reader->GetOutput() );
    histogramGenerator->SetNumberOfBins( numberOfBins );
    /*histogramGenerator->SetMarginalScale( 10.0 );
    histogramGenerator->SetHistogramMin( rangeCalculator->GetMinimum() - binWidth / 2 );
    histogramGenerator->SetHistogramMax( rangeCalculator->GetMaximum() + binWidth / 2 );*/
    histogramGenerator->Compute();
    typedef HistogramGeneratorType::HistogramType  HistogramType;
    HistogramType::ConstPointer histogram = histogramGenerator->GetOutput();

    // std::cout << "Histogram size " << histogram->Size() << std::endl;
    // for ( int binNumber = 0; binNumber < histogram->Size(); binNumber++ )
    //   {
    //   std::cout << "   binNumber " << binNumber << " -> frequency: "
    //             << histogram->GetFrequency( binNumber, 0 ) << std::endl;
    //   }

    const float  percentile = 0.05;
    const double  lowerBound = histogram->Quantile( 0, percentile );
    const double  upperBound = histogram->Quantile( 0, 1 - percentile );
    std::cout << "Lower " << percentile * 100 << "% percentile is at " << lowerBound << std::endl;
    std::cout << "Upper " << ( 1 - percentile ) * 100 << "% percentile is at " << upperBound << std::endl;


    // Reduce numbers of intensity levels. The filter works by setting a vector of
    // thresholds [ t0, t1, ..., tN ]; intensities below t0 map to 0, intensities
    // between t_i and t_{i+1} map to (i+1); and intensities over tN map to N+1
    std::vector< unsigned short >  thresholds;
    // const float  slope = static_cast< float >( rangeCalculator->GetMaximum() - rangeCalculator->GetMinimum() ) /
    //                       static_cast< float >( numberOfIntensityLevels );
    // const float  offset = rangeCalculator->GetMinimum();
    const float  slope = static_cast< float >( upperBound - lowerBound ) /
                         static_cast< float >( numberOfIntensityLevels );
    const float  offset = lowerBound;
    for ( int i = 1; i < numberOfIntensityLevels; i++ )
    {
      thresholds.push_back( static_cast< unsigned short >( slope * i + offset + 0.5 ) );
    }

    typedef itk::ThresholdLabelerImageFilter< ImageType, ImageType >  ReducerType;
    ReducerType::Pointer  reducer = ReducerType::New();
    reducer->SetThresholds( thresholds );
    reducer->SetInput( reader->GetOutput() );
    reducer->Update();

    // Write reduced image out
    typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer  writer = WriterType::New();
    writer->SetInput( reducer->GetOutput() );
    writer->SetFileName( "reduced.img" );
    writer->Write();


    // If we have tag-along images, apply the same transform to them as well
    if ( argc > 3 )
    {
      for ( int i = 3; i < argc; i++ )
      {
        // Read the image
        ReaderType::Pointer  reader = ReaderType::New();
        reader->SetFileName( argv[ i ] );
        reader->Update();

        // Reduce
        ReducerType::Pointer  reducer = ReducerType::New();
        reducer->SetThresholds( thresholds );
        reducer->SetInput( reader->GetOutput() );
        reducer->Update();

        // Write reduced image out
        std::ostringstream  fileNameStream;
        fileNameStream << "tagAlongReduced" << i-2 << ".img";
        WriterType::Pointer  writer = WriterType::New();
        writer->SetInput( reducer->GetOutput() );
        writer->SetFileName( fileNameStream.str().c_str() );
        writer->Write();
      }

    } // End test if there are tag-along images


  }
  catch( itk::ExceptionObject& e )
  {
    std::cerr << e << std::endl;
  }

  return 0;
};

