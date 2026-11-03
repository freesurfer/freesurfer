#ifndef __kvlAtlasMeshDeformationOptimizer_h
#define __kvlAtlasMeshDeformationOptimizer_h

#include "kvlAtlasMesh.h"
#include "kvlAtlasMeshPositionCostAndGradientCalculator.h"


namespace kvl
{

// Events generated
itkEventMacro( DeformationStartEvent, itk::UserEvent );
itkEventMacro( DeformationIterationEvent, itk::UserEvent );
itkEventMacro( DeformationEndEvent, itk::UserEvent );


/**
 *
 * Base class for (gradient-based) optimization of mesh deformation.
 * This and the various derived classes are based on the book
 * "Numerical Optimization" by Nocedal and Wright (Springer, 1999),
 * especially Chapter 3 for the line search implemented here,
 * Chapter 5 for the conjugate gradient subclass, and Chapter 9 for
 * the limited-memory BFGS subclass.
 * 
 */
class AtlasMeshDeformationOptimizer: public itk::Object
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshDeformationOptimizer  Self;
  typedef itk::Object  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  //itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshDeformationOptimizer, itk::Object );

  /** */
  void  SetMesh( AtlasMesh* mesh, 
                 bool resetIterations=true, 
                 bool resetMemory=true, 
                 bool reinitialize=true )
    {
    m_Mesh = mesh;
    if ( resetIterations )
      {
      m_IterationNumber = 0;
      }
    if ( resetMemory )
      {
      this->WipeMemory();
      }
    if ( reinitialize )
      {
      m_Initialized = false;
      }
      
    }

  /** */
  const AtlasMesh*  GetMesh() const
    {
    return m_Mesh;
    }

  //
  void  SetCostAndGradientCalculator( AtlasMeshPositionCostAndGradientCalculator* calculator, 
                                      bool resetIterations=false, 
                                      bool resetMemory=false, 
                                      bool reinitialize=false )
    {
    m_Calculator =  calculator;
    if ( resetIterations )
      {
      m_IterationNumber = 0;
      }
    if ( resetMemory )
      {
      this->WipeMemory();
      }
    if ( reinitialize )
      {
      m_Initialized = false;
      }
    }
    
  const AtlasMeshPositionCostAndGradientCalculator*  GetCostAndGradientCalculator() const
    {
    return m_Calculator;
    }
    
  // 
  unsigned int  GetIterationNumber() const
    { return m_IterationNumber; }
    
  //
  unsigned int  GetMaximumNumberOfIterations() const
    { return m_MaximumNumberOfIterations; }

  void  SetMaximumNumberOfIterations( unsigned int maximumNumberOfIterations )
    { m_MaximumNumberOfIterations = maximumNumberOfIterations; }

  //
  void  SetMaximalDeformationStopCriterion( double maximalDeformationStopCriterion )
    { m_MaximalDeformationStopCriterion = maximalDeformationStopCriterion; }

  //
  double  GetMaximalDeformationStopCriterion() const
    { return m_MaximalDeformationStopCriterion; }

  //
  unsigned int  GetIterationEventResolution() const
    { return m_IterationEventResolution; }
    
  //
  void SetIterationEventResolution( unsigned int  iterationEventResolution )
    { m_IterationEventResolution = iterationEventResolution; }
 
  /** */
  double GetMinLogLikelihoodTimesPrior() const
    {
    return m_Cost;  
    }
  
  /** */
  bool Go();
  
  /** */
  double Step();
  
  //
  void  SetVerbose( bool verbose )
    { m_Verbose = verbose; }
    
  //
  bool  GetVerbose() const
    { return m_Verbose; }

  //
  void  SetLineSearchMaximalDeformationIntervalStopCriterion( 
                           double lineSearchMaximalDeformationIntervalStopCriterion )
    {
    m_LineSearchMaximalDeformationIntervalStopCriterion = lineSearchMaximalDeformationIntervalStopCriterion;
    }
    
  //
  double  GetLineSearchMaximalDeformationIntervalStopCriterion() const
    {
    return m_LineSearchMaximalDeformationIntervalStopCriterion;
    }
    
    
protected:
  AtlasMeshDeformationOptimizer();
  virtual ~AtlasMeshDeformationOptimizer();

  /** */
  virtual double FindAndOptimizeNewSearchDirection() = 0;
 
  //
  virtual void Initialize();

  //
  virtual void WipeMemory()
    {      
    }

  //
  virtual void  GetCostAndGradient( const AtlasMesh::PointsContainer* position, 
                                    double& cost, 
                                    AtlasPositionGradientContainerType::Pointer& gradient );
  
  //
  double  ComputeMaximalDeformation( const AtlasPositionGradientContainerType* deformation ) const;
  
  // Compute position + alpha * deformationDirection
  void  AddDeformation( const AtlasMesh::PointsContainer* position, 
                        double alpha,
                        const AtlasPositionGradientContainerType* deformationDirection,                    
                        AtlasMesh::PointsContainer::Pointer&  newPosition,
                        double&  maximalDeformation ) const;
  
  // Compute inner product deformation1' * deformation2
  double  ComputeInnerProduct(  const AtlasPositionGradientContainerType* deformation1,
                                const AtlasPositionGradientContainerType* deformation2 ) const;
  
  // Compute beta1 * deformation1 + beta2 * deformation2
  AtlasPositionGradientContainerType::Pointer  
     LinearlyCombineDeformations( const AtlasPositionGradientContainerType* deformation1,
                                  double beta1,
                                  const AtlasPositionGradientContainerType* deformation2,
                                  double beta2 ) const;

  // Compute beta * deformation 
  AtlasPositionGradientContainerType::Pointer  
     ScaleDeformation( const AtlasPositionGradientContainerType* deformation,
                       double beta ) const;

  //
  void  DoLineSearch( const AtlasMesh::PointsContainer*  startPosition, 
                      double  startCost,
                      const AtlasPositionGradientContainerType*  startGradient,                    
                      const AtlasPositionGradientContainerType*  searchDirection,                    
                      double  startAlpha,
                      double  c1,
                      double  c2,
                      AtlasMesh::PointsContainer::Pointer&  newPosition,
                      double&  newCost,
                      AtlasPositionGradientContainerType::Pointer& newGradient,
                      double&  alphaUsed );

  //
  bool  m_Verbose;
  double  m_Cost;
  AtlasMesh::PointsContainer::Pointer  m_Position;
  AtlasPositionGradientContainerType::Pointer  m_Gradient;
  
  
private:
  AtlasMeshDeformationOptimizer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  
  //
  int  m_IterationNumber;
  int  m_MaximumNumberOfIterations;
  int  m_IterationEventResolution;
  bool  m_Initialized;


  AtlasMesh::Pointer  m_Mesh;
  AtlasMeshPositionCostAndGradientCalculator::Pointer  m_Calculator;
  double  m_MaximalDeformationStopCriterion;
  
  double  m_LineSearchMaximalDeformationLimit;
  double  m_LineSearchMaximalDeformationIntervalStopCriterion;
  
  
};


} // end namespace kvl

#endif

