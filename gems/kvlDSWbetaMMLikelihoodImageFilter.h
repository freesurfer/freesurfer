#ifndef kvlDSWbetaMMLikelihoodImageFilter_h
#define kvlDSWbetaMMLikelihoodImageFilter_h

#include "itkImageToImageFilter.h"
#include "kvlLikelihoodImageFilterBase.h"
#include "kvlGMMLikelihoodImageFilter.h"
#include "itkArray.h"
#include <vector>


namespace kvl
{
/**
 *
 * Input is four or more images. One or more of the first images may be intensity images
 * used for a Gaussian mixture model. These first intensity images are passed to a member
 * GMMlikelihoodImageFilter, which calculates likelihoods based on bunch of means and covariances.
 * The remaining seven images correspond to the FA and three elements of the leading eigenvector in
 * a DTI tensor image.
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
class DSWbetaMMLikelihoodImageFilter :
        public LikelihoodImageFilterBase< TInputImage >
{
public:

    typedef DSWbetaMMLikelihoodImageFilter    Self;
    typedef itk::SmartPointer< Self >        Pointer;
    typedef itk::SmartPointer< const Self >  ConstPointer;
    itkNewMacro(Self);

    itkTypeMacro(DSWbetaMMLikelihoodImageFilter, LikelihoodImageFilterBase);

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
                        const std::vector< double >& DSWbetaAlpha,
                        const std::vector< double >& DSWbetaBeta,
                        const std::vector< double >& DSWbetaVariance,
                        const std::vector< vnl_vector< double > >& DSWbetaMeans,
                        const std::vector< double >& negLogKummerSamples,
                        const double& negLogKummerIncrement,
                        const std::vector< double >&  DSWbetaMixtureWeights,
                        const std::vector< int >&  numberOfDSWbetaeniusPerClass,
                        const double& voxratio);

    //void SetGaussianImages( const std::vector<ImageType::ConstPointer>& images);

protected:
    DSWbetaMMLikelihoodImageFilter();

    virtual void BeforeThreadedGenerateData();

    virtual void ThreadedGenerateData(const RegionType & outputRegionForThread, itk::ThreadIdType);

private:
    DSWbetaMMLikelihoodImageFilter(const Self &);
    void operator=(const Self &);

    std::vector< double >  m_Alphas;
    std::vector< double >  m_Betas;
    std::vector< double >  m_LogBetaFunctions;
    std::vector< double >  m_Concentrations;
    std::vector< vnl_vector< double > >  m_means;
    std::vector< double >  m_logNormaliserSamples;
    double m_logNormaliserIncrement;
    std::vector< double >  m_DSWbetaMixtureWeights;
    std::vector< int >  m_NumberOfDSWbetaPerClass;
    double m_voxratio;

    //typedef GMMLikelihoodImageFilter< InputImageType >  GMMLikelihoodFilterType;
    //typename GMMLikelihoodFilterType::Pointer  m_likelihoodFilter;

};


}


#include "kvlDSWbetaMMLikelihoodImageFilter.hxx"

#endif // kvlDSWbetaMMLikelihoodImageFilter_h
