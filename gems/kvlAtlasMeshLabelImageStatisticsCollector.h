#ifndef __kvlAtlasMeshLabelImageStatisticsCollector_h
#define __kvlAtlasMeshLabelImageStatisticsCollector_h

#include "kvlAtlasMeshStatisticsCollector.h"
#include "kvlCompressionLookupTable.h"

namespace kvl
{


class AtlasMeshLabelImageStatisticsCollector: public AtlasMeshStatisticsCollector
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshLabelImageStatisticsCollector  Self;
  typedef AtlasMeshStatisticsCollector Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshLabelImageStatisticsCollector, itk::Object );

  /** Some typedefs */
  typedef CompressionLookupTable::ImageType  LabelImageType;

  /** */  
  void SetLabelImage( const LabelImageType* labelImage,
                      const CompressionLookupTable* lookupTable );

  
protected:
  AtlasMeshLabelImageStatisticsCollector();
  virtual ~AtlasMeshLabelImageStatisticsCollector();
  
  //  
  void GetContributionOfTetrahedron( const AtlasMesh::PointType& p0,
                                     const AtlasMesh::PointType& p1,
                                     const AtlasMesh::PointType& p2,
                                     const AtlasMesh::PointType& p3,
                                     const AtlasAlphasType&  alphasInVertex0,
                                     const AtlasAlphasType&  alphasInVertex1,
                                     const AtlasAlphasType&  alphasInVertex2,
                                     const AtlasAlphasType&  alphasInVertex3,
                                     double&  minLogLikelihood,
                                     AtlasAlphasType&  statisticsInVertex0,
                                     AtlasAlphasType&  statisticsInVertex1,
                                     AtlasAlphasType&  statisticsInVertex2,
                                     AtlasAlphasType&  statisticsInVertex3 );
   

private:
  AtlasMeshLabelImageStatisticsCollector(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  //
  LabelImageType::ConstPointer  m_LabelImage;
  CompressionLookupTable::ConstPointer  m_CompressionLookupTable;

};


} // end namespace kvl

#endif
