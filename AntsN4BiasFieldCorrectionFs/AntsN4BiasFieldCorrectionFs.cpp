/*
  This program wraps utilities developed within the ANTs toolbox.

        https://github.com/ANTsX/ANTs

  The itkN4BiasFieldCorrectionImageFilter class implements an nonuniform
  normalization algorithm described in the following paper.
  ITK-based spatially-adaptive filter described in the following paper.

        N. Tustison et al., N4ITK:  Improved N3 Bias Correction,
        IEEE Transactions on Medical Imaging, 29(6):1310-1320, June 2010.

  For ANTs license information, see distribution/docs/license.ants.txt
*/

#include "itkBSplineControlPointImageFilter.h"
#include "itkExpImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkN4BiasFieldCorrectionImageFilter.h"
#include "itkShrinkImageFilter.h"

#include "argparse.h"
#include "mri.h"
#include "AntsN4BiasFieldCorrectionFs.help.xml.h"
#ifdef _OPENMP
#include "romp_support.h"
#endif


int main(int argc, char **argv)
{
  // parse args
  ArgumentParser parser;
  parser.addHelp(AntsN4BiasFieldCorrectionFs_help_xml, AntsN4BiasFieldCorrectionFs_help_xml_len);
  parser.addArgument("-i", "--input",  1, String, true);
  parser.addArgument("-o", "--output", 1, String, true);
  parser.addArgument("-s", "--shrink", 1, Int, false);
  parser.addArgument("-t", "--iters", '+', Int, false);
  parser.addArgument("-d", "--dtype", 1, String, false);
  parser.addArgument("-r", "--replace-zeros", 3, Float, false);// offset scale remaskflag
  parser.addArgument("-e", "--threads-nondetermistic", 1, Int, false);
  parser.parse(argc, argv);

  std::string inputname = parser.retrieve<std::string>("input");
  std::string outputname = parser.retrieve<std::string>("output");
  std::vector<float> replace0 = parser.retrieve<std::vector<float>>("replace-zeros");

  // More than 1 thread is nondet in ITK
  int threads = 1;
  if(parser.exists("threads-nondetermistic")) threads = parser.retrieve<int>("threads-nondetermistic");
  #ifdef _OPENMP
  omp_set_num_threads(threads);
  #endif
  if(threads > 1) printf("threads %d > 1 -- ITK routines may be non-deterministic\n",threads);

  // read target data type (default is float)
  int dtype = MRI_FLOAT;
  if (parser.exists("dtype")) {
    std::string dtype_str = parser.retrieve<std::string>("dtype");
    std::transform(dtype_str.begin(), dtype_str.end(), dtype_str.begin(), ::tolower);
    if (dtype_str == "float") {
      dtype = MRI_FLOAT;
    } else if (dtype_str == "uchar") {
      dtype = MRI_UCHAR;
    } else if (dtype_str == "int") {
      dtype = MRI_INT;
    } else {
      fs::fatal() << "unrecognized target dtype '" << dtype_str << "'";
    }
  }

  MRI* mri = MRIread(inputname.c_str());

  MRI *masknonzero = NULL, *mritmp = NULL;
  int nreplace = 0;
  if(replace0.size()>0){
    // Replace 0s with small random values. This can be important when
    // the input is defaced.
    printf("Replacing zeros: offset=%f scale=%f\n",replace0[0],replace0[1]);
    if(mri->type != MRI_FLOAT){
      printf("Changing type from %d to MRI_FLOAT\n",mri->type);
      mritmp = MRIchangeType(mri,MRI_FLOAT,0,1,1);
      MRIfree(&mri);
      mri = mritmp;
    }
    masknonzero = MRIbinarize(mri, NULL, .000001, 0, 1);
    // Don't parallize this or won't be deterministic
    srand48(53); // seed for drand48()
    int c,r,s;
    for(c=0; c < mri->width; c++){
      for(r=0; r < mri->height; r++){
	for(s=0; s < mri->depth; s++){
	  if(MRIgetVoxVal(masknonzero,c,r,s,0)>0.5) continue;
	  double val = replace0[0] + replace0[1]*drand48();
	  MRIsetVoxVal(mri,c,r,s,0,val);
	  nreplace++;
	}
      }
    }
    printf("Replaced %d voxels\n",nreplace);
  }

  if (mri->nframes != 1) fs::fatal() << "input cannot be 4D (has " << mri->nframes << " frames)";

  // convert to ITK image
  ITKImageType::Pointer inputImage = mri->toITKImage();

  // set itk threads to 1
  itk::MultiThreader::SetGlobalDefaultNumberOfThreads(threads);

  // create a "mask" that just covers the whole image
  ITKImageType::Pointer maskImage = ITKImageType::New();
  maskImage->CopyInformation(inputImage);
  maskImage->SetRegions(inputImage->GetRequestedRegion());
  maskImage->Allocate(false);
  maskImage->FillBuffer(itk::NumericTraits<ITKImageType::PixelType>::OneValue());

  // init the bias field correcter
  typedef itk::N4BiasFieldCorrectionImageFilter<ITKImageType, ITKImageType, ITKImageType> CorrecterType;
  CorrecterType::Pointer correcter = CorrecterType::New();

  // convergence options
  // set number of iterations (default is 50x50x50x50)
  std::vector<int> numIters = {50, 50, 50, 50};
  if (parser.exists("iters")) numIters = parser.retrieve<std::vector<int>>("iters");
  CorrecterType::VariableSizeArrayType maximumNumberOfIterations(numIters.size());
  for (unsigned int d = 0; d < numIters.size(); d++) maximumNumberOfIterations[d] = numIters[d];
  correcter->SetMaximumNumberOfIterations(maximumNumberOfIterations);
  correcter->SetNumberOfFittingLevels(numIters.size());
  correcter->SetConvergenceThreshold(0.0);

  // shrink the image to save time
  int shrinkFactor = parser.exists("shrink") ? parser.retrieve<int>("shrink") : 4;

  typedef itk::ShrinkImageFilter<ITKImageType, ITKImageType> ShrinkerType;
  ShrinkerType::Pointer shrinker = ShrinkerType::New();
  shrinker->SetInput(inputImage);
  shrinker->SetShrinkFactors(shrinkFactor);

  // shrink the mask
  ShrinkerType::Pointer maskshrinker = ShrinkerType::New();
  maskshrinker->SetInput(maskImage);
  maskshrinker->SetShrinkFactors(shrinkFactor);

  // run the correction
  correcter->SetInput(shrinker->GetOutput());
  correcter->SetMaskImage(maskshrinker->GetOutput());
  correcter->Update();

  // reconstruct the bias field at full image resolution
  typedef itk::BSplineControlPointImageFilter<CorrecterType::BiasFieldControlPointLatticeType, CorrecterType::ScalarImageType> BSplinerType;
  BSplinerType::Pointer bspliner = BSplinerType::New();
  bspliner->SetInput(correcter->GetLogBiasFieldControlPointLattice());
  bspliner->SetSplineOrder(correcter->GetSplineOrder());
  bspliner->SetSize(inputImage->GetLargestPossibleRegion().GetSize());
  bspliner->SetOrigin(inputImage->GetOrigin());
  bspliner->SetDirection(inputImage->GetDirection());
  bspliner->SetSpacing(inputImage->GetSpacing());
  bspliner->Update();

  // extract log field from bspline
  ITKImageType::Pointer logField = ITKImageType::New();
  logField->SetLargestPossibleRegion(inputImage->GetLargestPossibleRegion());
  logField->SetBufferedRegion(inputImage->GetBufferedRegion());
  logField->SetRequestedRegion(inputImage->GetRequestedRegion());
  logField->SetSpacing(inputImage->GetSpacing());
  logField->SetOrigin(inputImage->GetOrigin());
  logField->SetDirection(inputImage->GetDirection());
  logField->Allocate();

  itk::ImageRegionIterator<CorrecterType::ScalarImageType> ItB(bspliner->GetOutput(), bspliner->GetOutput()->GetLargestPossibleRegion());
  itk::ImageRegionIterator<ITKImageType> ItF(logField, logField->GetLargestPossibleRegion());
  for(ItB.GoToBegin(), ItF.GoToBegin(); !ItB.IsAtEnd(); ++ItB, ++ItF) ItF.Set(ItB.Get()[0]);

  // exp bias field
  typedef itk::ExpImageFilter<ITKImageType, ITKImageType> ExpFilterType;
  ExpFilterType::Pointer expFilter = ExpFilterType::New();
  expFilter->SetInput(logField);
  expFilter->Update();

  // divide original image by the bias field to get corrected image
  typedef itk::DivideImageFilter<ITKImageType, ITKImageType, ITKImageType> DividerType;
  DividerType::Pointer divider = DividerType::New();
  divider->SetInput1(inputImage);
  divider->SetInput2(expFilter->GetOutput());
  divider->Update();

  // crop divider to the input size
  ITKImageType::RegionType inputRegion;
  inputRegion.SetIndex(inputImage->GetLargestPossibleRegion().GetIndex());
  inputRegion.SetSize(inputImage->GetLargestPossibleRegion().GetSize());

  typedef itk::ExtractImageFilter<ITKImageType, ITKImageType> CropperType;
  CropperType::Pointer cropper = CropperType::New();
  cropper->SetInput(divider->GetOutput());
  cropper->SetExtractionRegion(inputRegion);
  cropper->SetDirectionCollapseToSubmatrix();
  cropper->Update();

  // load ITK image back into MRI and write to disk
  MRI* dest = MRIallocSequence(mri->width, mri->height, mri->depth, dtype, mri->nframes);
  MRIcopyHeader(mri, dest);
  dest->loadITKImage(cropper->GetOutput());

  if(replace0.size()>0){
    if(nreplace > 0 && replace0[2]){
      printf("Remasking %d\n",nreplace);
      MRImask(dest, masknonzero, dest, 0, 0);
    }
    else {
      if(nreplace > 0) printf("Not remasking\n");
      else printf("No voxels to remask\n");
    }
  }

  MRIwrite(dest, outputname.c_str());
  MRIfree(&dest);
  MRIfree(&mri);

  printf("AntsN4BiasFieldCorrectionFs done\n");

  return 0;
}
