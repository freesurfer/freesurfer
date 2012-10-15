/**
 * @file  kvlEvaluateQuantitatively.cxx
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



int main( int argc, char* argv[] )
{

  // Sanity check on input
  if ( argc < 4 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " manualSegmentation automatedSegmentation compressionLookupTable.txt [ overrideHippocampusLabel ignoreEmptySlices ] " << std::endl;
    exit( -1 );
  }

  // Collect the input arguments
  const std::string  manualSegmentationFileName( argv[ 1 ] );
  const std::string  automatedSegmentationFileName( argv[ 2 ] );
  const std::string  compressionLookupTableFileName( argv[ 3 ] );
  bool  overrideHippocampusLabel = false;
  bool  ignoreEmptySlices = false;
  if ( argc > 4 )
  {
    std::istringstream  overrideHippocampusLabelStream( argv[ 4 ] );
    overrideHippocampusLabelStream >> overrideHippocampusLabel;
  }
  if ( argc > 5 )
  {
    std::istringstream  ignoreEmptySlicesStream( argv[ 5 ] );
    ignoreEmptySlicesStream >> ignoreEmptySlices;
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


    // Loop over all image voxels, collecting for each structure the manual volume, the automated volume, and
    // the overlap volume (needed for Dice coefficient)
    const int  numberOfLabels = compressionLookupTable->GetLabelStringLookupTable().size();
    std::cout << "Having " << numberOfLabels << " labels" << std::endl;
    itk::Array< float >  manualVolumes( numberOfLabels );
    manualVolumes.Fill( 0.0f );
    itk::Array< float >  automatedVolumes( numberOfLabels );
    automatedVolumes.Fill( 0.0f );
    itk::Array< float >  overlapVolumes( numberOfLabels );
    overlapVolumes.Fill( 0.0f );

    // Loop over all slices, ignoring those where nothing is know if the user so desires
    ImageType::IndexType  index = manualSegmentation->GetBufferedRegion().GetIndex();
    ImageType::SizeType  size = manualSegmentation->GetBufferedRegion().GetSize();
    const int  numberOfSlices = size[ 2 ];
    int  numberOfNonEmptySlices = 0;
    for ( int sliceNumber = 0; sliceNumber < numberOfSlices; sliceNumber++ )
    {
      // Construct region corresponding to this slice
      index[ 2 ] = sliceNumber;
      size[ 2 ] = 1;
      ImageType::RegionType  sliceRegion( index, size );

      //
      itk::Array< float >  sliceManualVolumes( numberOfLabels );
      sliceManualVolumes.Fill( 0.0f );
      itk::Array< float >  sliceAutomatedVolumes( numberOfLabels );
      sliceAutomatedVolumes.Fill( 0.0f );
      itk::Array< float >  sliceOverlapVolumes( numberOfLabels );
      sliceOverlapVolumes.Fill( 0.0f );


      itk::ImageRegionConstIterator< ImageType >  manualIt( manualSegmentation, sliceRegion );
      itk::ImageRegionConstIterator< ImageType >  automatedIt( automatedSegmentation, sliceRegion );
      for ( ; !manualIt.IsAtEnd(); ++manualIt, ++automatedIt )
      {
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

        sliceManualVolumes[ manualLabel ]++;
        sliceAutomatedVolumes[ automatedLabel ]++;

        if ( manualLabel == automatedLabel )
        {
          sliceOverlapVolumes[ manualLabel ]++;
        }

      } // End loop over all pixels in slice

      // Check if one of the slices is empty
      bool  manualIsEmpty = true;
      for ( unsigned int i = 1; i < sliceManualVolumes.size(); i++ )
      {
        if ( ( sliceManualVolumes[ i ] != 0 ) )
        {
          manualIsEmpty = false;
          break;
        }
      }
      bool  automatedIsEmpty = true;
      for ( unsigned int i = 1; i < sliceManualVolumes.size(); i++ )
      {
        if ( ( sliceAutomatedVolumes[ i ] != 0 ) )
        {
          automatedIsEmpty = false;
          break;
        }
      }

      const bool  isEmpty = ( manualIsEmpty || automatedIsEmpty );


      if ( ( !isEmpty ) || (!ignoreEmptySlices) )
      {
        if ( !isEmpty )
        {
          numberOfNonEmptySlices++;
        }

        // Add contributions of this slice to global statistics
        for ( unsigned int i = 0; i < sliceManualVolumes.size(); i++ )
        {
          manualVolumes[ i ] += sliceManualVolumes[ i ];
          automatedVolumes[ i ] += sliceAutomatedVolumes[ i ];
          overlapVolumes[ i ] += sliceOverlapVolumes[ i ];
        }

      } // End test if we should just skip this slice


    } // End loop over all slices


    // Print out
    std::cout << "Format: labelNumber   label    manualVolume   automatedVolume   Dice " << std::endl;
    for ( int labelNumber = 0; labelNumber < numberOfLabels; ++labelNumber )
    {
      std::cout << labelNumber << " -> "
                << compressionLookupTable->GetLabelStringLookupTable().find( labelNumber )->second << "    "
                << manualVolumes[ labelNumber ] << "    "
                << automatedVolumes[ labelNumber ] << "    "
                << 2.0f * ( overlapVolumes[ labelNumber ] ) / ( manualVolumes[ labelNumber ] + automatedVolumes[ labelNumber ] )
                << std::endl;
    }
    if ( ignoreEmptySlices )
    {
      std::cout << "Ignored all slices except for " << numberOfNonEmptySlices << std::endl;
    }


  }
  catch( itk::ExceptionObject& e )
  {
    std::cerr << e << std::endl;
  }

  return 0;
};
