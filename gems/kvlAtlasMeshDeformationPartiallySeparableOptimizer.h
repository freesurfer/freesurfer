#ifndef __kvlAtlasMeshDeformationPartiallySeparableSOptimizer_h
#define __kvlAtlasMeshDeformationPartiallySeparableSOptimizer_h

#include "kvlAtlasMeshDeformationOptimizer.h"
#include "vnl/vnl_matrix_fixed.h"


namespace kvl
{



class AtlasMeshDeformationPartiallySeparableSOptimizer: public AtlasMeshDeformationOptimizer
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshDeformationPartiallySeparableSOptimizer  Self;
  typedef AtlasMeshDeformationOptimizer  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshDeformationPartiallySeparableSOptimizer, AtlasMeshDeformationOptimizer );

  
protected:
  AtlasMeshDeformationPartiallySeparableSOptimizer();
  virtual ~AtlasMeshDeformationPartiallySeparableSOptimizer();
  
  void Initialize();

  double FindAndOptimizeNewSearchDirection(); 
  
private:
  AtlasMeshDeformationPartiallySeparableSOptimizer(const Self&); //purposely not implemented
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
