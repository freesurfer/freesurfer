#ifndef kvlFrobMMLikelihoodImageFilter_h
#define kvlFrobMMLikelihoodImageFilter_h

#include "itkImageToImageFilter.h"
#include "kvlLikelihoodImageFilterBase.h"
#include "kvlGMMLikelihoodImageFilter.h"
#include "itkArray.h"
#include <vector>


namespace kvl
{
/**
 *
 * Input is seven or more images. One or more of the first images may be intensity images
 * used for a Gaussian mixture model. These first intensity images are passed to a member
 * GMMlikelihoodImageFilter, which calculates likelihoods based on bunch of means and covariances.
 * The remaining seven images correspond to the log determinant and six independent entries in
 * a DTI tensor image. Based on Wishart degrees of freedom and scale matrices, combined with
 * the output of the GMM filter, likelihoods are calculated for each class as output.
 *
 * Design is taken from itkComposeImageFilter. I'm forcing the output to be of the form
 * itk::Image< itk::Array<xxx>, yyy> instead of the more general outputs that itkComposeImageFilter
 * can handle (such as itkVectorImage), because I known in advance most of pixels will be background
 * which I can then encode by assigning an itkArray of size 0 there (saving both memory and the need
 * for some other signal that no information was present to compute the class likelihoods in specific
 * voxels)
 *
 */

template< typename TInputImage >
class FrobMMLikelihoodImageFilter :
        public LikelihoodImageFilterBase< TInputImage >
{
public:

    typedef FrobMMLikelihoodImageFilter    Self;
    typedef itk::SmartPointer< Self >        Pointer;
    typedef itk::SmartPointer< const Self >  ConstPointer;
    itkNewMacro(Self);

    itkTypeMacro(FrobMMLikelihoodImageFilter, LikelihoodImageFilterBase);

    itkStaticConstMacro(Dimension, unsigned int, TInputImage::ImageDimension);

    typedef TInputImage                          InputImageType;
    typedef typename InputImageType::PixelType   InputPixelType;
    typedef itk::Image< itk::Array< InputPixelType >, TInputImage::ImageDimension >   OutputImageType;
    typedef typename OutputImageType::PixelType  OutputPixelType;
    typedef typename itk::NumericTraits< OutputPixelType >::ValueType  OutputPixelValueType;
    typedef typename InputImageType::RegionType  RegionType;
    typedef LikelihoodImageFilterBase< InputImageType >  Superclass;


    /** Some typedefs */
    typedef itk::Image< float, 3 >  ImageType;

    /** */
    void SetParameters( const int numberOfContrasts,
                        const std::vector< double >& frobVariance,
                        const std::vector< vnl_vector< double > >& frobMeans,
                        const std::vector< double >&  frobMixtureWeights,
                        const std::vector< int >&  numberOfFrobeniusPerClass,
                        const double voxratio);

    //void SetGaussianImages( const std::vector<ImageType::ConstPointer>& images);

protected:
    FrobMMLikelihoodImageFilter();

    virtual void BeforeThreadedGenerateData();

    virtual void ThreadedGenerateData(const RegionType & outputRegionForThread, itk::ThreadIdType);

private:
    FrobMMLikelihoodImageFilter(const Self &);
    void operator=(const Self &);

    std::vector< double >  m_Precisions;
    std::vector< vnl_vector< double > >  m_means;
    std::vector< double >  m_logNormaliser;
    std::vector< double >  m_frobMixtureWeights;
    std::vector< int >  m_NumberOfFrobeniusPerClass;
    double m_voxratio;

    //typedef GMMLikelihoodImageFilter< InputImageType >  GMMLikelihoodFilterType;
    //typename GMMLikelihoodFilterType::Pointer  m_likelihoodFilter;

};


}


#include "kvlFrobMMLikelihoodImageFilter.hxx"

#endif // kvlFrobMMLikelihoodImageFilter_h
