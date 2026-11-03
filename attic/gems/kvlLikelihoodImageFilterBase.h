#ifndef kvlLikelihoodImageFilterBase_h
#define kvlLikelihoodImageFilterBase_h

#include "itkImageToImageFilter.h"
#include "itkArray.h"
#include <vector>

namespace kvl
{
/** 
 * 
 * Input is one or more intensity images. Based on bunch of means and covariances, computes
 * likelihoods for each class as output.
 * 
 * Design is taken from itkComposeImageFilter. I'm forcing the output to be of the form 
 * itk::Image< itk::Array<xxx>, yyy> instead of the more general outputs that itkComposeImageFilter
 * can handle (such as itkVectorImage), because I known in advance most of pixels will be background
 * which I can then encode by assigning an itkArray of size 0 there (saving both memory and the need
 * for some other signal that no information was present to compute the class likelihoods in specific
 * voxels)
 * 
 * See https://itk.org/Doxygen/html/classitk_1_1ImageToImageFilter.html for itk::ImageToImageFilter 
 * implementation, and explanation on BeforeThreadedGenerateData() and ThreadedGenerateData().
 */

template< typename TInputImage >
class LikelihoodImageFilterBase:
  public itk::ImageToImageFilter< TInputImage,
                                  itk::Image< itk::Array< typename TInputImage::PixelType >,
                                              TInputImage::ImageDimension > >
{
public:

  typedef LikelihoodImageFilterBase    Self;
  typedef itk::SmartPointer< Self >        Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;
  itkNewMacro(Self);

  itkTypeMacro(LikelihoodImageFilterBase, itk::ImageToImageFilter);

  itkStaticConstMacro(Dimension, unsigned int, TInputImage::ImageDimension);

  typedef TInputImage                          InputImageType;
  typedef typename InputImageType::PixelType   InputPixelType;
  typedef itk::Image< itk::Array< InputPixelType >, TInputImage::ImageDimension >   OutputImageType;
  typedef typename OutputImageType::PixelType  OutputPixelType;
  typedef typename itk::NumericTraits< OutputPixelType >::ValueType  OutputPixelValueType;
  typedef typename InputImageType::RegionType  RegionType;
  typedef itk::ImageToImageFilter< InputImageType, OutputImageType >  Superclass;

protected:
  LikelihoodImageFilterBase();
  virtual ~LikelihoodImageFilterBase();

private:
  LikelihoodImageFilterBase(const Self &);
  void operator=(const Self &);

};


}

#include "kvlLikelihoodImageFilterBase.hxx"

#endif // kvlLikelihoodImageFilterBase_h
