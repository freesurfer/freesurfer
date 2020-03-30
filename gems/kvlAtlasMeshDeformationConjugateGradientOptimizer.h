#ifndef __kvlAtlasMeshDeformationConjugateGradientOptimizer_h
#define __kvlAtlasMeshDeformationConjugateGradientOptimizer_h

#include "kvlAtlasMeshDeformationOptimizer.h"


namespace kvl
{



class AtlasMeshDeformationConjugateGradientOptimizer: public AtlasMeshDeformationOptimizer
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshDeformationConjugateGradientOptimizer  Self;
  typedef AtlasMeshDeformationOptimizer  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshDeformationConjugateGradientOptimizer, AtlasMeshDeformationOptimizer );

protected:
  AtlasMeshDeformationConjugateGradientOptimizer();
  virtual ~AtlasMeshDeformationConjugateGradientOptimizer();
  
  void Initialize();

  double FindAndOptimizeNewSearchDirection(); 
  
private:
  AtlasMeshDeformationConjugateGradientOptimizer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  double  m_OldCost;
  AtlasPositionGradientContainerType::Pointer  m_OldGradient;
  AtlasPositionGradientContainerType::Pointer  m_OldSearchDirection;
  double  m_AlphaUsedLastTime;
  
  double  m_StartDistance;
  
  
};


} // end namespace kvl

#endif
