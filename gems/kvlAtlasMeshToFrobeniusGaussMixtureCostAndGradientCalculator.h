#ifndef kvlAtlasMeshToFrobeniusGaussMixtureCostAndGradientCalculator_h
#define kvlAtlasMeshToFrobeniusGaussMixtureCostAndGradientCalculator_h

#include "kvlAtlasMeshToIntensityImageCostAndGradientCalculatorBase.h"
#include "itkImage.h"
#include "kvlGMMLikelihoodImageFilter.h"
#include "kvlFrobMMLikelihoodImageFilter.h"


namespace kvl
{


class AtlasMeshToFrobeniusGaussMixtureCostAndGradientCalculator :
        public AtlasMeshToIntensityImageCostAndGradientCalculatorBase
{
public:

    /** Standard class typedefs */
    typedef AtlasMeshToFrobeniusGaussMixtureCostAndGradientCalculator  Self;
    typedef AtlasMeshToIntensityImageCostAndGradientCalculatorBase Superclass;
    typedef itk::SmartPointer< Self >  Pointer;
    typedef itk::SmartPointer< const Self >  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro( Self );

    /** Run-time type information (and related methods). */
    itkTypeMacro( AtlasMeshToFrobeniusGaussMixtureCostAndGradientCalculator,
                  kvlAtlasMeshToIntensityImageCostAndGradientCalculatorBase );

    /** Some typedefs */
    typedef itk::Image< float, 3 >  ImageType;

    /** */
    void SetParameters( const std::vector< vnl_vector< double > >& means,
                        const std::vector< vnl_matrix< double > >& variances,
                        const std::vector< double >&  mixtureWeights,
                        const std::vector< int >&  numberOfGaussiansPerClass );

    /** */
    void SetDiffusionParameters( const int numberOfContrasts,
                                 const std::vector< double >&  frobMixtureWeights,
                                 const std::vector< int >&  numberOfFrobeniusPerClass,
                                 const double& voxratio,
                                 const std::vector< double >& frobVariance,
                                 const std::vector< vnl_vector< double > >& frobMeans ) override;

    //void SetGaussianImages( const std::vector<ImageType::ConstPointer>& images);

    //void SetFrobeniusImages( const std::vector<ImageType::ConstPointer>& images);

    /** */
    //void Rasterize( const AtlasMesh* mesh ) override;

protected:
    AtlasMeshToFrobeniusGaussMixtureCostAndGradientCalculator();
    virtual ~AtlasMeshToFrobeniusGaussMixtureCostAndGradientCalculator();

    void AddDataContributionOfTetrahedron( const AtlasMesh::PointType& p0,
                                               const AtlasMesh::PointType& p1,
                                               const AtlasMesh::PointType& p2,
                                               const AtlasMesh::PointType& p3,
                                               const AtlasAlphasType&  alphasInVertex0,
                                               const AtlasAlphasType&  alphasInVertex1,
                                               const AtlasAlphasType&  alphasInVertex2,
                                               const AtlasAlphasType&  alphasInVertex3,
                                               ThreadAccumDataType&  priorPlusDataCost,
                                               AtlasPositionGradientThreadAccumType&  gradientInVertex0,
                                               AtlasPositionGradientThreadAccumType&  gradientInVertex1,
                                               AtlasPositionGradientThreadAccumType&  gradientInVertex2,
                                               AtlasPositionGradientThreadAccumType&  gradientInVertex3 ) override;

private:
  AtlasMeshToFrobeniusGaussMixtureCostAndGradientCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  //
  typedef GMMLikelihoodImageFilter< ImageType >  gmmLikelihoodFilterType;
  typedef FrobMMLikelihoodImageFilter< ImageType >  frobmmLikelihoodFilterType;
  //frobmmLikelihoodFilterType::Pointer  m_FrobMMLikelihoodFilter;
};


} // end namespace kvl

#endif // kvlAtlasMeshToFrobeniusGaussMixtureCostAndGradientCalculator_h
