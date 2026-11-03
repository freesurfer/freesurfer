#ifndef __kvlAtlasMeshToIntensityImageCostAndGradientCalculator_h
#define __kvlAtlasMeshToIntensityImageCostAndGradientCalculator_h

#include "kvlAtlasMeshToIntensityImageCostAndGradientCalculatorBase.h"
#include "itkImage.h"
#include "kvlGMMLikelihoodImageFilter.h"


namespace kvl
{


class AtlasMeshToIntensityImageCostAndGradientCalculator:
        public AtlasMeshToIntensityImageCostAndGradientCalculatorBase
{
public :
  
#ifdef GEMS_DEBUG_RASTERIZE_VOXEL_COUNT
  int m_tetrahedronCnt;
  int m_zeroVoxel_tetrahedronCnt;

  int m_totalVoxel;
  int m_totalVoxelInTetrahedron;
#endif
    
  /** Standard class typedefs */
  typedef AtlasMeshToIntensityImageCostAndGradientCalculator  Self;
  typedef AtlasMeshToIntensityImageCostAndGradientCalculatorBase Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshToIntensityImageCostAndGradientCalculator, kvlAtlasMeshToIntensityImageCostAndGradientCalculatorBase );

  /** Some typedefs */
  typedef itk::Image< float, 3 >  ImageType;

  /** */
  void SetParameters( const std::vector< vnl_vector< double > >& means, 
                      const std::vector< vnl_matrix< double > >& variances,
                      const std::vector< double >&  mixtureWeights,
                      const std::vector< int >&  numberOfGaussiansPerClass );

  
protected:
  AtlasMeshToIntensityImageCostAndGradientCalculator();
  virtual ~AtlasMeshToIntensityImageCostAndGradientCalculator();

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
  AtlasMeshToIntensityImageCostAndGradientCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  //
  typedef GMMLikelihoodImageFilter< ImageType >  gmmLikelihoodFilterType;
  
};


} // end namespace kvl

#endif
