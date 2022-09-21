#ifndef kvlLikelihoodImageFilterBase_h
#define kvlLikelihoodImageFilterBase_h

#include "itkImageToImageFilter.h"
#include "itkArray.h"
#include <vector>

namespace kvl
{

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
