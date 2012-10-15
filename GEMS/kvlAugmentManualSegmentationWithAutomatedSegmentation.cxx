/**
 * @file  kvlAugmentManualSegmentationWithAutomatedSegmentation.cxx
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
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"



int main( int argc, char** argv )
{
  // Sanity check on the input
  if ( argc < 3 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " automaticFileName manualFileName" << std::endl;
    return -1;
  }

  // Collect input parameters
  const std::string  automaticFileName = argv[ 1 ];
  const std::string  manualFileName = argv[ 2 ];

  //
  try
  {
    // Read the files in
    typedef itk::Image< unsigned short, 3 >  ImageType;
    typedef itk::ImageFileReader< ImageType >  ReaderType;
    ReaderType::Pointer  automaticReader = ReaderType::New();
    automaticReader->SetFileName( automaticFileName.c_str() );
    automaticReader->Update();
    ReaderType::Pointer  manualReader = ReaderType::New();
    manualReader->SetFileName( manualFileName.c_str() );
    manualReader->Update();

    // Loop over all voxels, and change labels in the manual image under certain conditions
    itk::ImageRegionConstIterator< ImageType >  automaticIt( automaticReader->GetOutput(),
        automaticReader->GetOutput()->GetBufferedRegion() );
    itk::ImageRegionIterator< ImageType >  manualIt( manualReader->GetOutput(),
        manualReader->GetOutput()->GetBufferedRegion() );
    for ( ; !automaticIt.IsAtEnd(); ++automaticIt, ++manualIt )
    {
      //std::cout << "Pixel has manual label of " << manualIt.Value() << std::endl;
      switch ( manualIt.Value() )
      {
        //case 44:  // Right-Inf-Lat-Vent
      case 53:  // Right-Hippocampus
        //case 63:  // Right-choroid-plexus
      case 500: // right_CA2-3
      case 502: // right_CA1
      case 503: // right_fimbria
      case 504: // right_presubiculum
      case 505: // right_hippocampal_fissure
      case 506: // right_CA4-DG
      case 507: // right_subiculum
      case 41:  // Right-Cerebral-White-Matter
        // Don't do anything; just keep the label as it was
        //std::cout << "     Keeping it as it is" << std::endl;
        break;
      default:
        // Look at the automatic segmentation and give it the correct FreeSurfer label
        switch ( automaticIt.Value() )
        {
        case 1:
          // CSF
          //std::cout << "    Changing it to automated CSF" << std::endl;
          manualIt.Value() = 24;
          break;
        case 2:
          // GM
          //std::cout << "    Changing it to automated GM" << std::endl;
          manualIt.Value() = 42;
          break;
        case 3:
          // WM
          //std::cout << "    Changing it to automated WM" << std::endl;
          manualIt.Value() = 41;
          break;
        }
      }

    } // End loop over all voxels

    // Write corrected file out
    typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer  writer = WriterType::New();
    writer->SetInput( manualReader->GetOutput() );
    writer->SetFileName( "combined.img" );
    writer->Write();

  }
  catch( itk::ExceptionObject& e )
  {
    std::cerr << e << std::endl;
  }

  return 0;
};

