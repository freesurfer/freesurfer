#include <iostream>
#include <string>         
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "GetPot.h"
#include "itkThresholdMaximumConnectedComponentsImageFilter.h"
#include "itkCoherenceEnhancingDiffusionImageFilter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryImageToLabelMapFilter.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkLabelOverlayImageFilter.h"
#include "mris_multimodal_refinement.h"

#include "mri.h"
#include "mrisurf.h"
#include "fsenv.h"
#include "cma.h"


int main(int narg, char * arg[])
{
	GetPot cl(narg, const_cast<char**>(arg));
	if(cl.size()==1 || cl.search(2,"--help","-h"))
	{
		std::cout<<"Usage: " << std::endl;
		std::cout<< arg[0] << " -t1 t1.mgz -t2 t2.mgz -aseg aseg.mgz -o output.mgz --shape "  << std::endl;   
		return -1;
	}



	const char *imageNameT1= cl.follow ("", "-t1");
	const char *imageNameT2 = cl.follow ("", "-t2");
	const char *outputName = cl.follow ("", "-o");

	if(!cl.search("--shape"))
	{
		MRI* imageAseg = NULL;
		if(cl.search("-aseg"))
		{
			const char *imageNameAseg = cl.follow ("", "-aseg");
			imageAseg =  MRIread(imageNameAseg) ;
		}


		MRI *imageT1 =  MRIread(imageNameT1) ;
		MRI *imageT2 =  MRIread(imageNameT2) ;
		MRI *vesselMR =  MRIcopy(imageT1, NULL) ;
		//MRI *whiteMR=  MRIcopy(imageT1, NULL) ;
		MRIS_MultimodalRefinement refinement;
		refinement.SegmentVessel(imageT1, imageT2, vesselMR, 0);

		MRIwrite(vesselMR,outputName) ;
		MRIfree(&imageT1);	
		MRIfree(&imageT2);	
		MRIfree(&vesselMR);	
	}
	else
	{
		typedef unsigned short LabelType;
		typedef itk::ShapeLabelObject< long unsigned int, 3> ShapeLabelObjectType;
		typedef itk::LabelMap< ShapeLabelObjectType > LabelMapType;

		typedef itk::Image<unsigned char, 3> ImageType;
		typedef itk::Image<float, 3> FloatImageType;
	
		typedef itk::CastImageFilter<FloatImageType, ImageType > CastToCharType;

		itk::ImageFileReader<FloatImageType>::Pointer reader = itk::ImageFileReader<FloatImageType>::New();
		reader->SetFileName(imageNameT1);
		reader->Update();
		using BinaryFilterType = itk::BinaryThresholdImageFilter<FloatImageType, FloatImageType >;	
		BinaryFilterType::Pointer binarize = BinaryFilterType::New();  
		binarize->SetLowerThreshold( 0 );
		binarize->SetUpperThreshold(1);
		binarize->SetInput(reader->GetOutput());
		binarize->SetInsideValue(255);
		binarize->SetOutsideValue(0);
		binarize->Update();

		CastToCharType::Pointer castFilter2 = CastToCharType::New();
		castFilter2->SetInput(binarize->GetOutput());
		castFilter2->Update();	


		using BinaryImageToLabelMapFilterType = itk::BinaryImageToShapeLabelMapFilter<ImageType>;
		BinaryImageToLabelMapFilterType::Pointer binaryImageToLabelMapFilter = BinaryImageToLabelMapFilterType::New();
		binaryImageToLabelMapFilter->SetInput(castFilter2->GetOutput());
		binaryImageToLabelMapFilter->Update();
		// The output of this filter is an itk::LabelMap, which contains itk::LabelObject's
		std::cout << "There are " << binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() << " objects." << std::endl;

		LabelMapType *labelMap = binaryImageToLabelMapFilter->GetOutput();
		std::cout  <<  labelMap->GetNumberOfLabelObjects() << " labels." << std::endl;

		std::vector<unsigned long> labelsToRemove;

		// Retrieve all attributes
		for (unsigned int n = 0; n < labelMap->GetNumberOfLabelObjects(); ++n)
		{
			ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(n);
			if (labelObject->GetPhysicalSize() >1000)
			{
				labelsToRemove.push_back(labelObject->GetLabel());
			}
			else if (labelObject->GetRoundness() > .80)
			{
				if( labelObject->GetPhysicalSize() > 10)
				{
					labelObject->SetLabel(1);
				}
				else
				{
					labelObject->SetLabel(2);
				}
			}
			else if (labelObject->GetElongation() > .99)
			{

				labelObject->SetLabel(3);
			}
		}

		// Remove all regions that were marked for removal.
		for(unsigned int i = 0; i < labelsToRemove.size(); ++i)
		{
			binaryImageToLabelMapFilter->GetOutput()->RemoveLabel(labelsToRemove[i]);
		}

		typedef itk::LabelMapToLabelImageFilter<BinaryImageToLabelMapFilterType::OutputImageType, ImageType> LabelMapToLabelImageFilterType;
		LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
		labelMapToLabelImageFilter->SetInput(binaryImageToLabelMapFilter->GetOutput());
		labelMapToLabelImageFilter->Update();

		typedef itk::RGBPixel<unsigned char> RGBPixelType;
		typedef itk::Image<RGBPixelType, 3> RGBImageType;

		typedef itk::LabelOverlayImageFilter<FloatImageType, ImageType, RGBImageType> LabelOverlayImageFilterType;
		LabelOverlayImageFilterType::Pointer labelOverlayImageFilter = LabelOverlayImageFilterType::New();
		labelOverlayImageFilter->SetInput(reader->GetOutput());
		labelOverlayImageFilter->SetLabelImage(labelMapToLabelImageFilter->GetOutput());
		labelOverlayImageFilter->SetOpacity(.5);
		labelOverlayImageFilter->Update();

		typedef  itk::ImageFileWriter< RGBImageType  > WriterTypeS;
		WriterTypeS::Pointer writers = WriterTypeS::New();
		writers->SetFileName(outputName);
		writers->SetInput(labelOverlayImageFilter->GetOutput());
		writers->Update();


	}
	return 0;
}
