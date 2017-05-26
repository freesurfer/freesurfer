#ifndef __kvlMutualInformationCostAndGradientCalculator_h
#define __kvlMutualInformationCostAndGradientCalculator_h

#include "kvlAtlasMeshPositionCostAndGradientCalculator.h"
#include "kvlHistogrammer.h"


namespace kvl
{


/**
 *
 */
class MutualInformationCostAndGradientCalculator: public AtlasMeshPositionCostAndGradientCalculator
{
public :
  
  /** Standard class typedefs */
  typedef MutualInformationCostAndGradientCalculator  Self;
  typedef AtlasMeshPositionCostAndGradientCalculator Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( MutualInformationCostAndGradientCalculator, AtlasMeshPositionCostAndGradientCalculator );

  /** Some typedefs */
  typedef itk::Image< float, 3 >  ImageType;

  /** */  
  void SetImage( const ImageType* image );

  /** */  
  void Rasterize( const AtlasMesh* mesh );
  
  
protected:
  MutualInformationCostAndGradientCalculator();
  virtual ~MutualInformationCostAndGradientCalculator();
  
  void AddDataContributionOfTetrahedron( const AtlasMesh::PointType& p0,
                                         const AtlasMesh::PointType& p1,
                                         const AtlasMesh::PointType& p2,
                                         const AtlasMesh::PointType& p3,
                                         const AtlasAlphasType&  alphasInVertex0,
                                         const AtlasAlphasType&  alphasInVertex1,
                                         const AtlasAlphasType&  alphasInVertex2,
                                         const AtlasAlphasType&  alphasInVertex3,
                                         double&  priorPlusDataCost,
                                         AtlasPositionGradientType&  gradientInVertex0,
                                         AtlasPositionGradientType&  gradientInVertex1,
                                         AtlasPositionGradientType&  gradientInVertex2,
                                         AtlasPositionGradientType&  gradientInVertex3 );
  
private:
  MutualInformationCostAndGradientCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  //
  Histogrammer::Pointer  m_Histogrammer;
  
};


} // end namespace kvl

#endif
