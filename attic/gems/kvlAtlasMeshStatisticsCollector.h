#ifndef __kvlAtlasMeshStatisticsCollector_h
#define __kvlAtlasMeshStatisticsCollector_h

#include "kvlAtlasMeshRasterizor.h"

namespace kvl
{


class AtlasMeshStatisticsCollector: public AtlasMeshRasterizor
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshStatisticsCollector  Self;
  typedef AtlasMeshRasterizor Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshStatisticsCollector, itk::Object );

  /** Some typedefs */
  typedef itk::MapContainer< AtlasMesh::PointIdentifier , AtlasAlphasType >   StatisticsContainerType;

  /** */
  const StatisticsContainerType*  GetLabelStatistics() const
    {
    return m_LabelStatistics;
    }

  /** */
  double GetMinLogLikelihood() const
    {
    return m_MinLogLikelihood;
    }

  /** */  
  void Rasterize( const AtlasMesh* mesh );
  
  
protected:
  AtlasMeshStatisticsCollector();
  virtual ~AtlasMeshStatisticsCollector();
  
  //
  bool RasterizeTetrahedron( const AtlasMesh* mesh, 
                             AtlasMesh::CellIdentifier tetrahedronId,
                             int threadNumber );

  virtual void GetContributionOfTetrahedron( const AtlasMesh::PointType& p0,
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
                                             AtlasAlphasType&  statisticsInVertex3 )
   {
   }  

  
private:
  AtlasMeshStatisticsCollector(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  //
  StatisticsContainerType::Pointer  m_LabelStatistics;
  double  m_MinLogLikelihood;

   //
  std::vector< StatisticsContainerType::Pointer >  m_ThreadSpecificLabelStatistics;
  std::vector< double >  m_ThreadSpecificMinLogLikelihoods;

};


} // end namespace kvl

#endif
