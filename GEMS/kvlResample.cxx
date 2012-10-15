/**
 * @file  kvlResample.cxx
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
#include "itkRescaleIntensityImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkMinimumMaximumImageCalculator.h"


/*
* Note that when higher-order interpolators are used, the intensities in the resampled image can
* fall outside the intensity range of the original image, and indeed, outside the intensity range
* of the pixel type. Therefore, the resampling is done first with floating point precision, and the
* resulting image is then rescaled to fall within the intensity range of the original image.
*/




int main( int argc, char* argv[] )
{
  // Sanity check on input
  if ( argc < 3 )
  {
    std::cerr << argv[0] << " inputImage referenceImage [ interpolationMode=1 computeDefaultPixelValue=1 defaultPixelValue=0.0 ]" << std::endl;
    return -1;
  }

  // Parse input
  const std::string  inputImageFileName = argv[ 1 ];
  const std::string  referenceImageFileName = argv[ 2 ];

  unsigned int  interpolationMode = 1;
  if ( argc > 3 )
  {
    std::istringstream  interpolationModeStream( argv[ 3 ] );
    interpolationModeStream >> interpolationMode;
  }

  bool  computeDefaultPixelValue = true;
  if ( argc > 4 )
  {
    std::istringstream  computeDefaultPixelValueStream( argv[ 4 ] );
    computeDefaultPixelValueStream >> computeDefaultPixelValue;
  }

  float  defaultPixelValue = 0.0f;
  if ( argc > 5 )
  {
    std::istringstream  defaultPixelValueStream( argv[ 5 ] );
    defaultPixelValueStream >> defaultPixelValue;
  }


  // Add support for MGH file format to ITK. An alternative way to add this by default would be
  // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
  itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );

  try
  {
    // Some typedefs first...
    typedef short  PixelType;
    typedef itk::Image< PixelType, 3 >  ImageType;
    typedef float  InternalPixelType;
    typedef itk::Image< InternalPixelType, 3 >  InternalImageType;

    typedef itk::ImageFileReader< ImageType >  ReaderType;

    typedef itk::ResampleImageFilter< ImageType, InternalImageType >  ResamplerType;
    typedef ResamplerType::InterpolatorType  InterpolatorType;
    typedef ResamplerType::PixelType  DefaultPixelValueType;

    typedef itk::RescaleIntensityImageFilter< InternalImageType, ImageType >  RescalerType;
    typedef itk::CastImageFilter< InternalImageType, ImageType >   CasterType;

    typedef itk::NearestNeighborInterpolateImageFunction< ImageType, double > NNInterpolatorType;
    typedef itk::LinearInterpolateImageFunction< ImageType, double >   LinearInterpolatorType;
    typedef itk::BSplineInterpolateImageFunction< ImageType, double > BSplineInterpolatorType;

    typedef itk::MinimumMaximumImageCalculator< ImageType >  MinimumMaximumCalculatorType;

    typedef itk::ImageFileWriter< ImageType >  WriterType;


    // Read images
    ReaderType::Pointer  inputReader = ReaderType::New();
    inputReader->SetFileName( inputImageFileName.c_str() );
    inputReader->Update();
    ImageType::ConstPointer  inputImage = inputReader->GetOutput();

    ReaderType::Pointer  referenceReader = ReaderType::New();
    referenceReader->SetFileName( referenceImageFileName.c_str() );
    referenceReader->Update();
    ImageType::ConstPointer  referenceImage = referenceReader->GetOutput();


    // If desired, set the default pixel value to the minimum of the input image
    MinimumMaximumCalculatorType::Pointer  minimumMaximumCalculator = MinimumMaximumCalculatorType::New();
    minimumMaximumCalculator->SetImage( inputImage );
    if ( computeDefaultPixelValue )
    {
      minimumMaximumCalculator->Compute();
      defaultPixelValue =
        static_cast< DefaultPixelValueType >( minimumMaximumCalculator->GetMinimum() );

      std::cout << "Setting default pixel value to the minimum of the input image (calculated as "
                << defaultPixelValue<< ")" << std::endl;
    }


    // Set up the correct interpolator
    InterpolatorType::Pointer interpolator = 0;
    BSplineInterpolatorType::Pointer bsplineInterpolator = 0;
    switch( interpolationMode )
    {
    case 0 :
      std::cout << "Using Nearest Neighbor interpolation" << std::endl;
      interpolator = NNInterpolatorType::New();
      break;
    case 1 :
      std::cout << "Using Linear interpolation" << std::endl;
      interpolator = LinearInterpolatorType::New();
      break;
    default :
      std::cout << "Using Cubic B-Spline interpolation" << std::endl;
      bsplineInterpolator = BSplineInterpolatorType::New();
      bsplineInterpolator->SetSplineOrder( 3 );
      interpolator = bsplineInterpolator;
    }


    // Set up the resampler
    ResamplerType::Pointer resampler = ResamplerType::New();
    resampler->SetOutputParametersFromImage( referenceImage );
    resampler->SetInterpolator( interpolator );
    resampler->SetDefaultPixelValue( defaultPixelValue );
    resampler->SetInput( inputImage );
    kvl::ProgressReporter  reporter( resampler, "Resampling" );

    ImageType::Pointer  resampledImage = 0;
    if ( interpolationMode == 0 )
    {
      // Set up the caster
      CasterType::Pointer caster = CasterType::New();
      caster->SetInput( resampler->GetOutput() );
      kvl::ProgressReporter  reporter( caster, "Casting image back" );

      // Let the beast go!
      caster->Update();

      //
      resampledImage = caster->GetOutput();
    }
    else
    {
      //
      minimumMaximumCalculator->Compute();

      // Set up the rescaler
      RescalerType::Pointer rescaler = RescalerType::New();
      rescaler->SetInput( resampler->GetOutput() );
      rescaler->SetOutputMinimum( minimumMaximumCalculator->GetMinimum() /* itk::NumericTraits<PixelType>::min() */ );
      rescaler->SetOutputMaximum( minimumMaximumCalculator->GetMaximum() /* itk::NumericTraits<PixelType>::max() */ );
      kvl::ProgressReporter  reporter( rescaler, "Rescaling image" );

      // Let the beast go!
      rescaler->Update();

      //
      resampledImage = rescaler->GetOutput();
    }


    // Write out
    std::ostringstream  outputFileNameStream;
    outputFileNameStream << itksys::SystemTools::GetFilenameWithoutExtension( inputImageFileName.c_str() )
                         << "_resampled.mgz";
    const std::string  outputFileName = outputFileNameStream.str();
    WriterType::Pointer  writer = WriterType::New();
    writer->SetInput( resampledImage );
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

