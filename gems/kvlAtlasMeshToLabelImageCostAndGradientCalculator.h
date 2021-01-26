#ifndef __kvlAtlasMeshToLabelImageCostAndGradientCalculator_h
#define __kvlAtlasMeshToLabelImageCostAndGradientCalculator_h

#include "kvlAtlasMeshPositionCostAndGradientCalculator.h"
#include "kvlCompressionLookupTable.h"

namespace kvl
{


class AtlasMeshToLabelImageCostAndGradientCalculator: public AtlasMeshPositionCostAndGradientCalculator
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshToLabelImageCostAndGradientCalculator  Self;
  typedef AtlasMeshPositionCostAndGradientCalculator Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshToLabelImageCostAndGradientCalculator, AtlasMeshPositionCostAndGradientCalculator );

  // Some typedefs
  typedef CompressionLookupTable::ImageType  LabelImageType;
  
  // Set label image
  void SetLabelImage( const LabelImageType* labelImage, 
                      const CompressionLookupTable*  lookupTable );
  
protected:
  AtlasMeshToLabelImageCostAndGradientCalculator();
  virtual ~AtlasMeshToLabelImageCostAndGradientCalculator();
  
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
                                         AtlasPositionGradientThreadAccumType&  gradientInVertex3 );
  
private:
  AtlasMeshToLabelImageCostAndGradientCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  //
  LabelImageType::ConstPointer  m_LabelImage;
  CompressionLookupTable::ConstPointer  m_CompressionLookupTable;
  
};


} // end namespace kvl

#endif
