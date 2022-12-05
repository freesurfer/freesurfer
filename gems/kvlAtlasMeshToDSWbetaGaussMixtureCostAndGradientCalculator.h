#ifndef kvlAtlasMeshToDSWbetaGaussMixtureCostAndGradientCalculator_h
#define kvlAtlasMeshToDSWbetaGaussMixtureCostAndGradientCalculator_h

#include "kvlAtlasMeshToIntensityImageCostAndGradientCalculatorBase.h"
#include "itkImage.h"
#include "kvlGMMLikelihoodImageFilter.h"
#include "kvlDSWbetaMMLikelihoodImageFilter.h"


namespace kvl
{


class AtlasMeshToDSWbetaGaussMixtureCostAndGradientCalculator :
        public AtlasMeshToIntensityImageCostAndGradientCalculatorBase
{
public:

    /** Standard class typedefs */
    typedef AtlasMeshToDSWbetaGaussMixtureCostAndGradientCalculator  Self;
    typedef AtlasMeshToIntensityImageCostAndGradientCalculatorBase Superclass;
    typedef itk::SmartPointer< Self >  Pointer;
    typedef itk::SmartPointer< const Self >  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro( Self );

    /** Run-time type information (and related methods). */
    itkTypeMacro( AtlasMeshToDSWbetaGaussMixtureCostAndGradientCalculator,
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
                                 const std::vector< double >&  DSWbetaMixtureWeights,
                                 const std::vector< int >&  numberOfDSWbetaePerClass,
                                 const double& voxratio,
                                 const std::vector< double >& DSWbetaAlpha,
                                 const std::vector< vnl_vector< double > >& DSWbetaMeans,
                                 const std::vector< double >& DSWbetaBeta,
                                 const std::vector< double >& DSWbetaConcentration,
                                 const std::vector< double >& logKummerSamples,
                                 const double& logKummerIncrement ) override;

    //void SetGaussianImages( const std::vector<ImageType::ConstPointer>& images);

    //void SetDSWbetaImages( const std::vector<ImageType::ConstPointer>& images);

    /** */
    //void Rasterize( const AtlasMesh* mesh ) override;

protected:
    AtlasMeshToDSWbetaGaussMixtureCostAndGradientCalculator();
    virtual ~AtlasMeshToDSWbetaGaussMixtureCostAndGradientCalculator();

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
  AtlasMeshToDSWbetaGaussMixtureCostAndGradientCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  //
  typedef GMMLikelihoodImageFilter< ImageType >  gmmLikelihoodFilterType;
  typedef DSWbetaMMLikelihoodImageFilter< ImageType >  dswbetammLikelihoodFilterType;
  //dswbetammLikelihoodFilterType::Pointer  m_DSWbetaMMLikelihoodFilter;
};


} // end namespace kvl

#endif // kvlAtlasMeshToDSWbetaGaussMixtureCostAndGradientCalculator_h
