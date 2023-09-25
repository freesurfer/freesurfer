/*
  This program wraps utilities developed within the ANTs toolbox.

        https://github.com/ANTsX/ANTs

  The itkAdaptiveNonLocalMeansDenoisingImageFilter class (from ANTs) implements an
  ITK-based spatially-adaptive filter described in the following paper.

        J. V. Manjon, P. Coupe, Luis Marti-Bonmati, D. L. Collins, and M. Robles. Adaptive 
        Non-Local Means Denoising of MR Images With Spatially Varying Noise Levels, 
        Journal of Magnetic Resonance Imaging, 31:192-203, June 2010. 

  A technical description of the code is available at the following link.

        https://www.insight-journal.org/browse/publication/979

  For ANTs license information, see distribution/docs/license.ants.txt
*/


#include "itkAdaptiveNonLocalMeansDenoisingImageFilter.h"
#include "argparse.h"
#include "mri.h"
#include "AntsDenoiseImageFs.help.xml.h"


int main(int argc, char **argv)
{
  // parse args
  ArgumentParser parser;
  parser.addHelp(AntsDenoiseImageFs_help_xml, AntsDenoiseImageFs_help_xml_len);
  parser.addArgument("-i", "--input",  1, String, true);
  parser.addArgument("-o", "--output", 1, String, true);
  parser.addArgument("--rician");
  parser.parse(argc, argv);

  std::string inputname = parser.retrieve<std::string>("input");
  std::string outputname = parser.retrieve<std::string>("output");

  MRI* mri = MRIread(inputname.c_str());

  // set itk threads to 1
#if ITK_VERSION_MAJOR >= 5
  itk::MultiThreaderBase::SetGlobalDefaultNumberOfThreads(1);
#else
  itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);  
#endif  

  for (int frame = 0 ; frame < mri->nframes ; frame++) {
    // convert frame to ITK image
    ITKImageType::Pointer image = mri->toITKImage(frame);

    // configure denoiser
    typedef itk::AdaptiveNonLocalMeansDenoisingImageFilter<ITKImageType, ITKImageType> DenoiserType;
    typename DenoiserType::Pointer denoiser = DenoiserType::New();
    denoiser->SetInput(image);

    // set noise model
    bool rician = parser.exists("rician");
    denoiser->SetUseRicianNoiseModel(rician);

    // set neighborhood patch radius
    typename DenoiserType::NeighborhoodRadiusType neighborhoodPatchRadius;
    neighborhoodPatchRadius.Fill(1);
    denoiser->SetNeighborhoodPatchRadius(neighborhoodPatchRadius);

    // set neighborhood search radius
    typename DenoiserType::NeighborhoodRadiusType neighborhoodSearchRadius;
    neighborhoodSearchRadius.Fill(2);
    denoiser->SetNeighborhoodSearchRadius(neighborhoodSearchRadius);

    // set local mean and var neighborhood radius
    typename DenoiserType::NeighborhoodRadiusType neighborhoodRadiusForLocalMeanAndVariance;
    neighborhoodRadiusForLocalMeanAndVariance.Fill(1);
    denoiser->SetNeighborhoodRadiusForLocalMeanAndVariance(neighborhoodRadiusForLocalMeanAndVariance);

    // hardcoded options from the ANTs source
    denoiser->SetEpsilon(0.00001);
    denoiser->SetMeanThreshold(0.95);
    denoiser->SetVarianceThreshold(0.5);
    denoiser->SetSmoothingFactor(1.0);
    denoiser->SetSmoothingVariance(2.0);

    denoiser->Update();

    // run denoising and copy into MRI structure
    typename ITKImageType::Pointer denoised = denoiser->GetOutput();
    denoised->Update();
    denoised->DisconnectPipeline();
    mri->loadITKImage(denoised, frame);
  }

  MRIwrite(mri, outputname.c_str());
  MRIfree(&mri);

  return 0;
}

