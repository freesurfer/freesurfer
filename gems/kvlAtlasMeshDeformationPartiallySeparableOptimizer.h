#ifndef __kvlAtlasMeshDeformationPartiallySeparableOptimizer_h
#define __kvlAtlasMeshDeformationPartiallySeparableOptimizer_h

#include "kvlAtlasMeshDeformationOptimizer.h"
#include "vnl/vnl_matrix_fixed.h"


namespace kvl
{



class AtlasMeshDeformationPartiallySeparableOptimizer: public AtlasMeshDeformationOptimizer
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshDeformationPartiallySeparableOptimizer  Self;
  typedef AtlasMeshDeformationOptimizer  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshDeformationPartiallySeparableOptimizer, AtlasMeshDeformationOptimizer );

  
protected:
  AtlasMeshDeformationPartiallySeparableOptimizer();
  virtual ~AtlasMeshDeformationPartiallySeparableOptimizer();
  
  void Initialize();

  double FindAndOptimizeNewSearchDirection(); 
  
private:
  AtlasMeshDeformationPartiallySeparableOptimizer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  double  m_OldCost;
  AtlasPositionGradientContainerType::Pointer  m_OldGradient;
  AtlasPositionGradientContainerType::Pointer  m_OldSearchDirection;
  double  m_AlphaUsedLastTime;
  
  double  m_StartDistance;
  
  typedef vnl_matrix_fixed< double, 12, 12 >  miniApproxHessianType;
  std::vector< miniApproxHessianType >  m_MiniApproxHessians;
  
};


} // end namespace kvl

#endif
