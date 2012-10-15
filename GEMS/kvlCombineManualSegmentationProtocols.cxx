/**
 * @file  kvlCombineManualSegmentationProtocols.cxx
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
    std::cerr << "Usage: " << argv[ 0 ] << " fileNameOfProtocol1 fileNameOfProtocol2" << std::endl;
    return -1;
  }

  // Collect input parameters
  const std::string  fileNameOfProtocol1 = argv[ 1 ];
  const std::string  fileNameOfProtocol2 = argv[ 2 ];

  //
  try
  {
    // Read the files in
    typedef itk::Image< unsigned short, 3 >  ImageType;
    typedef itk::ImageFileReader< ImageType >  ReaderType;
    ReaderType::Pointer  reader1 = ReaderType::New();
    reader1->SetFileName( fileNameOfProtocol1.c_str() );
    reader1->Update();
    ReaderType::Pointer  reader2 = ReaderType::New();
    reader2->SetFileName( fileNameOfProtocol2.c_str() );
    reader2->Update();

    // Loop over all voxels, and change labels in image 2 under certain conditions
    itk::ImageRegionConstIterator< ImageType >  it1( reader1->GetOutput(),
        reader1->GetOutput()->GetBufferedRegion() );
    itk::ImageRegionIterator< ImageType >  it2( reader2->GetOutput(),
        reader2->GetOutput()->GetBufferedRegion() );
    for ( ; !it1.IsAtEnd(); ++it1, ++it2 )
    {
      if ( it2.Value() == 500 ) // right_CA2/3
      {
        if ( it1.Value() == 501 ) // right_alveus
        {
          it2.Value() = 501;
        }


      }

    }

    // Write corrected file out
    typedef itk::ImageFileWriter< ImageType >  WriterType;
    WriterType::Pointer  writer = WriterType::New();
    writer->SetInput( reader2->GetOutput() );
    writer->SetFileName( "combinedProtocol.mhd" );
    writer->Write();

  }
  catch( itk::ExceptionObject& e )
  {
    std::cerr << e << std::endl;
  }

  return 0;
};

