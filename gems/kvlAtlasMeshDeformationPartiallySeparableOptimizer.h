#ifndef __kvlAtlasMeshDeformationPartiallySeparableOptimizer_h
#define __kvlAtlasMeshDeformationPartiallySeparableOptimizer_h

#include "kvlAtlasMeshDeformationOptimizer.h"
#include "vnl/vnl_matrix_fixed.h"

#define KVL_ENABLE_MULTITREADING

#ifdef KVL_ENABLE_MULTITREADING
  #include "Eigen/Sparse"
#endif



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

  /** */
  void SetNumberOfThreads( int numberOfThreads )
    { 
    m_NumberOfThreads = numberOfThreads;
    }
  
  /** */
  int GetNumberOfThreads() const
    {
    return m_NumberOfThreads;
    }
  
  
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
  
  std::map< AtlasMesh::PointIdentifier, int >   m_NodeNumberLookupTable;
  std::map< AtlasMesh::CellIdentifier, int >  m_TetrahedronNumberLookupTable;
  
#ifdef KVL_ENABLE_MULTITREADING
  //
  int  m_NumberOfThreads;
  typedef Eigen::SparseMatrix< double >  HessianType;
  std::vector< HessianType >  m_ThreadSpecificHessians;
  static ITK_THREAD_RETURN_TYPE ThreaderCallback( void *arg );
  struct ThreadStruct
    {
    Pointer  m_Optimizer;
    AtlasPositionGradientContainerType::ConstPointer  m_S;
    AtlasPositionGradientContainerType::ConstPointer  m_Y;
    AtlasMesh::ConstPointer  m_Mesh;
    std::vector< AtlasMesh::CellIdentifier >  m_TetrahedronIds;
    };
#endif    

};


} // end namespace kvl

#endif
