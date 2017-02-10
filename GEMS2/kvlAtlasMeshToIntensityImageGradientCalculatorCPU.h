#ifndef __kvlAtlasMeshToIntensityImageGradientCalculatorCPU_h
#define __kvlAtlasMeshToIntensityImageGradientCalculatorCPU_h

#include "kvlAtlasMeshRasterizorCPU.h"
#include "itkImage.h"
#include "itkAffineTransform.h"
#include "kvlGaussianLikelihoodImageFilter.h"


namespace kvl
{


/**
 *
 */
class AtlasMeshToIntensityImageGradientCalculatorCPU: public AtlasMeshRasterizorCPU
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshToIntensityImageGradientCalculatorCPU  Self;
  typedef AtlasMeshRasterizorCPU Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshToIntensityImageGradientCalculatorCPU, itk::Object );

  /** Some typedefs */
  typedef itk::Image< float, 3 >  ImageType;
  typedef itk::AffineTransform< double, 3 >  TransformType;

  /** */  
  void SetImages( const std::vector< ImageType::Pointer >& images );

  /** */
  void SetMeshToImageTransform( const TransformType* meshToImageTransform );

  /** */
  void SetIgnoreDeformationPrior( bool ignoreDeformationPrior )
    { m_IgnoreDeformationPrior = ignoreDeformationPrior; }

  /** */
  void SetOnlyDeformationPrior( bool onlyDeformationPrior )
    { m_OnlyDeformationPrior = onlyDeformationPrior; }

  /** */
  const AtlasPositionGradientContainerType* GetPositionGradient() const
    {
    return m_PositionGradient;
    }

  /** */
  AtlasPositionGradientContainerType* GetPositionGradient()
    {
    return m_PositionGradient;
    }
  
  /** */
  double GetMinLogLikelihoodTimesPrior() const
    {
    return m_MinLogLikelihoodTimesPrior;
    }

  /** */
  void SetMeans( const std::vector< vnl_vector<float> >& means );

  //  
  void SetPrecisions( const std::vector< vnl_matrix<float> >& precisions );
    
  /** */  
  void Rasterize( const AtlasMesh* mesh );
  
  
protected:
  AtlasMeshToIntensityImageGradientCalculatorCPU();
  virtual ~AtlasMeshToIntensityImageGradientCalculatorCPU();
  
  //
  bool RasterizeTetrahedron( const AtlasMesh* mesh, 
                             AtlasMesh::CellIdentifier tetrahedronId,
                             int threadNumber );

private:
  AtlasMeshToIntensityImageGradientCalculatorCPU(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  //
  void ImposeSlidingBoundaryConditions( const AtlasMesh* mesh );
 
  //
  typedef GaussianLikelihoodImageFilter< ImageType >  LikelihoodFilterType;
  LikelihoodFilterType::Pointer  m_LikelihoodFilter;
  
  //
  AtlasPositionGradientContainerType::Pointer  m_PositionGradient;
  bool  m_IgnoreDeformationPrior;
  bool  m_OnlyDeformationPrior;
  double  m_MinLogLikelihoodTimesPrior;

  //
  typedef itk::Matrix< double >  SlidingBoundaryCorrectionMatrixType;
  SlidingBoundaryCorrectionMatrixType  m_SlidingBoundaryCorrectionMatrices[ 8 ]; 
  
  //
  std::vector< AtlasPositionGradientContainerType::Pointer >  m_ThreadSpecificPositionGradients;
  std::vector< double >  m_ThreadSpecificMinLogLikelihoodTimesPriors;
  bool  m_Abort;

};


} // end namespace kvl

#endif
