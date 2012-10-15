/**
 * @file  kvlEvaluateQuantitativelyWithHausdorffDistance.cxx
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
#include "kvlCompressionLookupTable.h"
#include "itkImageFileReader.h"
#include "itkImageRegionConstIterator.h"
#include "itkHausdorffDistanceImageFilter.h"
#include "itkContourDirectedMeanDistanceImageFilter.h"
#include "itkContourMeanDistanceImageFilter.h"


int main( int argc, char* argv[] )
{

  // Sanity check on input
  if ( argc < 4 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " manualSegmentation automatedSegmentation compressionLookupTable.txt [ overrideHippocampusLabel ] " << std::endl;
    exit( -1 );
  }

  // Collect the input arguments
  const std::string  manualSegmentationFileName( argv[ 1 ] );
  const std::string  automatedSegmentationFileName( argv[ 2 ] );
  const std::string  compressionLookupTableFileName( argv[ 3 ] );
  bool  overrideHippocampusLabel = false;
  if ( argc > 4 )
  {
    std::istringstream  overrideHippocampusLabelStream( argv[ 4 ] );
    overrideHippocampusLabelStream >> overrideHippocampusLabel;
  }


  try
  {
    // Read the manual segmentation
    typedef kvl::CompressionLookupTable::CompressedImageType  ImageType;
    typedef itk::ImageFileReader< ImageType >  ReaderType;
    ReaderType::Pointer  reader = ReaderType::New();
    reader->SetFileName( manualSegmentationFileName.c_str() );
    reader->Update();
    ImageType::Pointer  tmp = reader->GetOutput();
    tmp->DisconnectPipeline();
    ImageType::ConstPointer  manualSegmentation = tmp.GetPointer();

    // Read the automated segmentation
    reader->SetFileName( automatedSegmentationFileName.c_str() );
    reader->Update();
    ImageType::ConstPointer  automatedSegmentation = reader->GetOutput();

    // Read the compressionLookupTable
    kvl::CompressionLookupTable::Pointer  compressionLookupTable = kvl::CompressionLookupTable::New();
    if ( !compressionLookupTable->Read( compressionLookupTableFileName.c_str() ) )
    {
      std::cerr << "Couldn't read compressionLookupTable from file " << compressionLookupTableFileName << std::endl;
      exit( -1 );
    }


    // Look up to what label the FreeSurfer "53 Right-Hippocampus" maps, since we're possibly not gonna counts this
    // as an error because it's not reproducible across subjects in the manual delineations
    ImageType::PixelType  rightHippocampusLabel( 0 );
    kvl::CompressionLookupTable::CompressionLookupTableType::const_iterator  it = compressionLookupTable->GetCompressionLookupTable().find( 53 );
    if ( it != compressionLookupTable->GetCompressionLookupTable().end() )
    {
      rightHippocampusLabel = it->second;
      std::cout << "Found that Right-Hippocampus maps to labelNumber " << static_cast< int >( rightHippocampusLabel ) << std::endl;
    }
    else
    {
      std::cout << "Didn't find Right-Hippocampus in the lookup table " << std::endl;
      overrideHippocampusLabel = false;
    }



    // Loop over all labels
    std::cout << "Having " <<  compressionLookupTable->GetLabelStringLookupTable().size() << " labels" << std::endl;
    std::cout << "Format: labelNumber   label    manualVolume   automatedVolume   hausdorffDistance   directedMeanDistance   meanDistance" << std::endl;
    for ( kvl::CompressionLookupTable::LabelStringLookupTableType::const_iterator  it = compressionLookupTable->GetLabelStringLookupTable().begin();
          it != compressionLookupTable->GetLabelStringLookupTable().end(); ++it )
    {
      const int  labelNumber = it->first;
      const std::string  labelName = it->second;
      //std::cout << "Analyzing " << labelName << " which has labelNumber " << labelNumber << std::endl;

      // Create empty image for binarized manual and automated segmentations
      ImageType::Pointer  binaryManualSegmentation = ImageType::New();
      binaryManualSegmentation->SetRegions( manualSegmentation->GetBufferedRegion() );
      binaryManualSegmentation->Allocate();

      ImageType::Pointer  binaryAutomatedSegmentation = ImageType::New();
      binaryAutomatedSegmentation->SetRegions( manualSegmentation->GetBufferedRegion() );
      binaryAutomatedSegmentation->Allocate();

      // Fill in contents. Also count manual and automated volume as a double-check
      int  manualVolume = 0;
      int  automatedVolume = 0;
      itk::ImageRegionConstIterator< ImageType >  manualIt( manualSegmentation, manualSegmentation->GetBufferedRegion() );
      itk::ImageRegionConstIterator< ImageType >  automatedIt( automatedSegmentation, automatedSegmentation->GetBufferedRegion() );
      itk::ImageRegionIterator< ImageType >  binaryManualIt( binaryManualSegmentation, manualSegmentation->GetBufferedRegion() );
      itk::ImageRegionIterator< ImageType >  binaryAutomatedIt( binaryAutomatedSegmentation, automatedSegmentation->GetBufferedRegion() );
      for ( ; !manualIt.IsAtEnd(); ++manualIt, ++automatedIt, ++binaryManualIt, ++binaryAutomatedIt )
      {
        // Retrieve original labels and correct if necessary
        ImageType::PixelType  manualLabel =  manualIt.Value();
        ImageType::PixelType  automatedLabel =  automatedIt.Value();

        if ( overrideHippocampusLabel )
        {
          if ( manualLabel == rightHippocampusLabel )
          {
            automatedLabel = rightHippocampusLabel;
          }
          if ( automatedLabel == rightHippocampusLabel )
          {
            manualLabel = rightHippocampusLabel;
          }
        }

        // Fill in binary images
        if ( manualLabel == labelNumber )
        {
          binaryManualIt.Value() = 1;
          manualVolume++;
        }
        else
        {
          binaryManualIt.Value() = 0;
        }

        if ( automatedLabel == labelNumber )
        {
          binaryAutomatedIt.Value() = 1;
          automatedVolume++;
        }
        else
        {
          binaryAutomatedIt.Value() = 0;
        }


      } // End loop over all pixels

      // Now compute Hausdorff distance
      typedef itk::HausdorffDistanceImageFilter< ImageType, ImageType >  HausdorffCalculatorType;
      HausdorffCalculatorType::Pointer  hausdorffCalculator = HausdorffCalculatorType::New();
      hausdorffCalculator->SetInput1( binaryManualSegmentation );
      hausdorffCalculator->SetInput2( binaryAutomatedSegmentation );
      hausdorffCalculator->Update();
      const double  hausdorffDistance = hausdorffCalculator->GetHausdorffDistance();

      // Compute directed mean distance: loop over all boundary points in manual segmentation, compute
      // distance to automated surface, and average
      typedef itk::ContourDirectedMeanDistanceImageFilter< ImageType, ImageType > DirectedMeanCalculatorType;
      DirectedMeanCalculatorType::Pointer  directedMeanCalculator = DirectedMeanCalculatorType::New();
      directedMeanCalculator->SetInput1( binaryManualSegmentation );
      directedMeanCalculator->SetInput2( binaryAutomatedSegmentation );
      directedMeanCalculator->Update();
      const double  directedMeanDistance = directedMeanCalculator->GetContourDirectedMeanDistance();

      // Compute mean distance: calculate directed mean both ways, and retain maximum value
      typedef itk::ContourMeanDistanceImageFilter< ImageType, ImageType >  MeanCalculatorType;
      MeanCalculatorType::Pointer  meanCalculator = MeanCalculatorType::New();
      meanCalculator->SetInput1( binaryManualSegmentation );
      meanCalculator->SetInput2( binaryAutomatedSegmentation );
      meanCalculator->Update();
      const double  meanDistance = meanCalculator->GetMeanDistance();


      // Print out result
      std::cout << labelNumber << " -> "
                << labelName << "    "
                << manualVolume << "    "
                << automatedVolume << "    "
                << hausdorffDistance << "    "
                << directedMeanDistance << "    "
                << meanDistance << "    "
                << std::endl;

    } // End loop over all labels


  }
  catch( itk::ExceptionObject& e )
  {
    std::cerr << e << std::endl;
  }

  return 0;
};
