#ifndef __kvlAtlasMeshToIntensityImageLogDomainCostAndGradientCalculator_h
#define __kvlAtlasMeshToIntensityImageLogDomainCostAndGradientCalculator_h

#include "kvlAtlasMeshToIntensityImageCostAndGradientCalculator.h"
#include "itkImage.h"
#include "kvlGMMLikelihoodImageFilter.h"


namespace kvl
{


class AtlasMeshToIntensityImageLogDomainCostAndGradientCalculator:
        public AtlasMeshToIntensityImageCostAndGradientCalculator
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshToIntensityImageLogDomainCostAndGradientCalculator  Self;
  typedef AtlasMeshToIntensityImageCostAndGradientCalculator Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshToIntensityImageLogDomainCostAndGradientCalculator, kvlAtlasMeshToIntensityImageCostAndGradientCalculator );

  /** Some typedefs */
  typedef itk::Image< float, 3 >  ImageType;

  /** */
  //void SetParameters( const std::vector< vnl_vector< double > >& means, 
  //                    const std::vector< vnl_matrix< double > >& variances,
  //                    const std::vector< double >&  mixtureWeights,
  //                    const std::vector< int >&  numberOfGaussiansPerClass );

  
protected:
  AtlasMeshToIntensityImageLogDomainCostAndGradientCalculator();
  virtual ~AtlasMeshToIntensityImageLogDomainCostAndGradientCalculator();

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
  AtlasMeshToIntensityImageLogDomainCostAndGradientCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  //
  typedef GMMLikelihoodImageFilter< ImageType >  gmmLikelihoodFilterType;

  
};


} // end namespace kvl

#endif
