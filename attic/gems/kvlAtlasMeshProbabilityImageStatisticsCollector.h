#ifndef __kvlAtlasMeshProbabilityImageStatisticsCollector_h
#define __kvlAtlasMeshProbabilityImageStatisticsCollector_h

#include "kvlAtlasMeshStatisticsCollector.h"
#include "kvlCompressionLookupTable.h"

namespace kvl
{


class AtlasMeshProbabilityImageStatisticsCollector: public AtlasMeshStatisticsCollector
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshProbabilityImageStatisticsCollector  Self;
  typedef AtlasMeshStatisticsCollector Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshProbabilityImageStatisticsCollector, itk::Object );

  /** Some typedefs */
  typedef itk::Image< AtlasAlphasType, 3 >  ProbabilityImageType;

  /** */  
  void SetProbabilityImage( const ProbabilityImageType*  probabilityImage )
    {
    m_ProbabilityImage = probabilityImage;  
    }  

  
protected:
  AtlasMeshProbabilityImageStatisticsCollector();
  virtual ~AtlasMeshProbabilityImageStatisticsCollector();
  
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
  AtlasMeshProbabilityImageStatisticsCollector(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  //
  ProbabilityImageType::ConstPointer  m_ProbabilityImage;

};


} // end namespace kvl

#endif
