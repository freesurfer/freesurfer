#ifndef __kvlAtlasMeshDeformationFixedStepGradientDescentOptimizer_h
#define __kvlAtlasMeshDeformationFixedStepGradientDescentOptimizer_h

#include "kvlAtlasMeshDeformationOptimizer.h"


namespace kvl
{



class AtlasMeshDeformationFixedStepGradientDescentOptimizer: public AtlasMeshDeformationOptimizer
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshDeformationFixedStepGradientDescentOptimizer  Self;
  typedef AtlasMeshDeformationOptimizer  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshDeformationFixedStepGradientDescentOptimizer, AtlasMeshDeformationOptimizer );

  //
  void SetStepSize( double stepSize )
    {
    m_StepSize = stepSize;
    }  
  
  double GetStepSize() const
    {
    return m_StepSize;
    }
  

protected:
  AtlasMeshDeformationFixedStepGradientDescentOptimizer();
  virtual ~AtlasMeshDeformationFixedStepGradientDescentOptimizer();
  
  double FindAndOptimizeNewSearchDirection();
  
private:
  AtlasMeshDeformationFixedStepGradientDescentOptimizer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  double  m_StepSize;
  double  m_LineSearchStopCriterion;
  
};
 

} // end namespace kvl

#endif
