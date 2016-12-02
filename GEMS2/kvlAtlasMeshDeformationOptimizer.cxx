#include "kvlAtlasMeshDeformationOptimizer.h"


namespace kvl
{


//
//
//
AtlasMeshDeformationOptimizer
::AtlasMeshDeformationOptimizer()
{
  //m_Image = 0;
  m_ProbabilityImage = 0;
  m_SegmentedImage = 0;
  m_Mesh = 0;
  m_MeshToImageTransform = 0;

  m_Initialized = false;

  m_IterationNumber = 0;
  m_MaximumNumberOfIterations = 1000 /* itk::NumericTraits< unsigned int >::max() */;
  m_IterationEventResolution = 10;
  m_NumberOfThreads = itk::MultiThreader::GetGlobalDefaultNumberOfThreads();
  m_Verbose = true;

  m_mapCompToComp = 0;
}


//
//
//
AtlasMeshDeformationOptimizer
::~AtlasMeshDeformationOptimizer()
{
}


} // end namespace kvl

