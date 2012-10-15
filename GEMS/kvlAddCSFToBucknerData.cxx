/**
 * @file  kvlAddCSFToBucknerData.cxx
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
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMGHImageIOFactory.h"


int main( int argc, char* argv[] )
{

  // Sanity check on input
  if ( argc < 3 )
  {
    std::cerr << argv[0] << " originalSegmentation normImage" << std::endl;
    return -1;
  }

  // Add support for MGH file format to ITK. An alternative way to add this by default would be
  // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
  itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );

  // Read the segmentation from file
  typedef itk::Image< unsigned char, 3 >  ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  ReaderType::Pointer  reader = ReaderType::New();
  reader->SetFileName( argv[ 1 ] );
  reader->Update();
  ImageType::Pointer  segmentation = reader->GetOutput();
  segmentation->DisconnectPipeline();

  // Read the image from file
  reader->SetFileName( argv[ 2 ] );
  reader->Update();
  ImageType::ConstPointer  image = reader->GetOutput();


  // Loop over all voxels, and if a voxel has label 0 in the segmentation but
  // a non-zero intensity, replace the label by label 24 (FreeSurfer speak for CSF)
  itk::ImageRegionConstIterator< ImageType >  imageIt( image,
      image->GetBufferedRegion() );
  itk::ImageRegionIterator< ImageType >  segmentationIt( segmentation,
      segmentation->GetBufferedRegion() );
  for ( ; !imageIt.IsAtEnd(); ++imageIt, ++segmentationIt )
  {
    if ( ( segmentationIt.Value() == 0 ) && ( imageIt.Value() > 0 ) )
    {
      segmentationIt.Value() = 24; // CSF
    }
  }


  // Write result out
  const std::string  outputName = "BucknerSegmentationWithAddedCSF.mgz";
  typedef itk::ImageFileWriter< ImageType >  WriterType;
  WriterType::Pointer  writer = WriterType::New();
  writer->SetFileName( outputName.c_str() );
  writer->SetInput( segmentation );
  writer->Update();

  std::cout << "Wrote to file " << outputName << std::endl;


  return 0;
};

