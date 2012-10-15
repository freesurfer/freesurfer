/**
 * @file  kvlUpsample.cxx
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
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMGHImageIOFactory.h"
#include "kvlProgressReporter.h"
#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"


/*
 * Let's implement a "smart" image upsampler. Each voxel's intensity will be repeated several times for
 * upsamplingFactorX * upsamplingFactorY * upsamplingFactorZ subvoxels, but the origin is set so that
 * the center position of the original voxels are the *average* position of the corresponding subvoxels.
 *
 * This useful when modeling partial volume effects etc
 *
 */


/*
 * For the implementation, recall that the upsampledImageGridCoordinates "y" (i.e., position counted in
 * voxel indices) relates to the originalImageGridCoordinaes "x" as follows:
 *
 *     y = f * x + ( f - 1 ) / 2
 *
 * where f is the upsampling factor.
 *
 *
 * In matrix notation: y = diag( fx, fy, fz ) * x + ( (fx-1)/2 (fy-1)/2 (fz-1)/2 )^T
 *
 * On the other hand, the image-to-world transformation of the original image is given by
 *
 *    x_w = R * diag( s_x, s_y, s_z ) * x + ( o_x o_y o_z )^T
 *
 * Therefore, we have
 *
 *    x_w = R * diag( s_x, s_y, s_z ) * diag( 1/fx, 1/fy, 1/fz ) * ( y - ( (fx-1)/2 (fy-1)/2 (fz-1)/2 )^T ) + ( o_x o_y o_z )^T
 *
 * which yields
 *
 *    x_w = R * S_N * y + o_N
 *
 * as the resulting image-to-world transformation of the upsampled image,
 *
 * where
 *
 *     S_N = diag( s_x / fx, s_y / fy, s_z / fz )   contains the smaller voxel spacings
 *
 * and
 *
 *     o_N = ( o_x o_y o_z )^T  - R * ( s_x/fx*(fx-1)/2  s_y/fy*(fy-1)/2  s_z/fz*(fz-1)/2 )^T
 *
 */



/*
 * One more thing that is useful for the implementation is to know that the
 * first upsampledImageGridCoordinates "y" (i.e., position counted in
 * voxel indices) corresponding to originalImageGridCoordinaes "x" is given by:
 *
 *     z = f * ( x - 0.5 ) + ( f - 1 ) / 2 + 1/2
 *       = f * ( x - 1/2 ) + f / 2
 *       = f * x
 *
 * where f is the upsampling factor.
 *
 */


int main( int argc, char* argv[] )
{
  // Sanity check on input
  if ( argc < 3 )
  {
    std::cerr << argv[0] << " inputImage upsamplingFactorX [ upsamplingFactorY=upsamplingFactorX upsamplingFactorZ=upsamplingFactorX ]" << std::endl;
    return -1;
  }

  // Parse input
  const std::string  inputImageFileName = argv[ 1 ];
  int  upsamplingFactorX = 1;
  std::istringstream  upsamplingFactorXStream( argv[ 2 ] );
  upsamplingFactorXStream >> upsamplingFactorX;

  int  upsamplingFactorY = upsamplingFactorX;
  if ( argc > 3 )
  {
    std::istringstream  upsamplingFactorYStream( argv[ 3 ] );
    upsamplingFactorYStream >> upsamplingFactorY;
  }

  int  upsamplingFactorZ = upsamplingFactorX;
  if ( argc > 4 )
  {
    std::istringstream  upsamplingFactorZStream( argv[ 4 ] );
    upsamplingFactorZStream >> upsamplingFactorZ;
  }

  std::cout << "inputImage: " << inputImageFileName << std::endl;
  std::cout << "upsamplingFactorX: " << upsamplingFactorX << std::endl;
  std::cout << "upsamplingFactorY: " << upsamplingFactorY << std::endl;
  std::cout << "upsamplingFactorZ: " << upsamplingFactorZ << std::endl;
  const int upsamplingFactors[] = { upsamplingFactorX, upsamplingFactorY, upsamplingFactorZ };

  // Add support for MGH file format to ITK. An alternative way to add this by default would be
  // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
  itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );

  try
  {
    // Some typedefs first...
    typedef short  PixelType;
    typedef itk::Image< PixelType, 3 >  ImageType;

    typedef itk::ImageFileReader< ImageType >  ReaderType;

    typedef itk::ResampleImageFilter< ImageType, ImageType >  ResamplerType;

    typedef itk::NearestNeighborInterpolateImageFunction< ImageType, double > InterpolatorType;

    typedef itk::ImageFileWriter< ImageType >  WriterType;


    // Read image
    ReaderType::Pointer  reader = ReaderType::New();
    reader->SetFileName( inputImageFileName.c_str() );
    reader->Update();
    ImageType::ConstPointer  inputImage = reader->GetOutput();


    // Set up the resampler
    ResamplerType::Pointer resampler = ResamplerType::New();
    resampler->SetInput( inputImage );
    resampler->SetInterpolator( InterpolatorType::New() );
    ImageType::SpacingType  outputSpacing;
    ImageType::PointType  outputOrigin;
    ImageType::SizeType   outputSize;
    ImageType::IndexType  outputStartIndex;
    for ( int i = 0; i < 3; i++ )
    {
      outputSpacing[ i ] = inputImage->GetSpacing()[ i ] / static_cast< double >( upsamplingFactors[ i ] );
      outputOrigin[ i ] = outputSpacing[ i ] * ( static_cast< double >( upsamplingFactors[ i ] ) - 1.0 ) / 2.0;
      outputSize[ i ] = inputImage->GetLargestPossibleRegion().GetSize()[ i ] * upsamplingFactors[ i ];
      outputStartIndex[ i ] = inputImage->GetLargestPossibleRegion().GetIndex()[ i ] *  upsamplingFactors[ i ];
    }
    outputOrigin = inputImage->GetDirection() * outputOrigin;
    for ( int i = 0; i < 3; i++ )
    {
      outputOrigin[ i ] = inputImage->GetOrigin()[ i ] - outputOrigin[ i ];
    }
    resampler->SetOutputOrigin( outputOrigin );
    resampler->SetOutputSpacing( outputSpacing );
    resampler->SetOutputDirection( inputImage->GetDirection() );
    resampler->SetOutputStartIndex( outputStartIndex );
    resampler->SetSize( outputSize );
    kvl::ProgressReporter  reporter( resampler, "Resampling" );


    // Write out
    std::ostringstream  outputFileNameStream;
    outputFileNameStream << itksys::SystemTools::GetFilenameWithoutExtension( inputImageFileName.c_str() )
                         << "_upsampled_"
                         << upsamplingFactorX << "_"
                         << upsamplingFactorY << "_"
                         << upsamplingFactorZ << ".mgz";
    const std::string  outputFileName = outputFileNameStream.str();
    WriterType::Pointer  writer = WriterType::New();
    writer->SetInput( resampler->GetOutput() );
    writer->SetFileName( outputFileName.c_str() );
    writer->Update();

    std::cout << "Wrote to file " << outputFileName << std::endl;
  }
  catch ( itk::ExceptionObject& e )
  {
    std::cerr << e << std::endl;
    return -1;
  }


  return 0;
};

