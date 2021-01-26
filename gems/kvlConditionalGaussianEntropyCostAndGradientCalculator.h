#ifndef __kvlConditionalGaussianEntropyCostAndGradientCalculator_h
#define __kvlConditionalGaussianEntropyCostAndGradientCalculator_h

#include "kvlAtlasMeshPositionCostAndGradientCalculator.h"
#include "itkImage.h"


namespace kvl
{


/**
 *
 * Cost function closely related to Maes's Mutual Information (MI) and Roche's
 * correlation ratio: essentially computes conditional entropy H(X|Y) where Y
 * is a mesh-based probabilistic atlas and X is an intensity image, assuming
 * that each of the classes in the atlas have an associated Gaussian conditional
 * intensity distribution with unknown mean and variance, which are estimated 
 * non-iteratively for each mesh deformation in Maes's "partial volume 
 * interpolation" (PVI) fashion. (for comparison, Maes's MI is H(X) - H(X|Y) and
 * uses a categorial distribution on the intensities instead, but is otherwise
 * identical).
 * 
 * As per Maes's idea, the PVI trick allows for a background class to be ignored
 * in the cost function computation, with voxels fractionally appearing and 
 * disappearing at the image/mesh borders as the mesh moves/deforms; this guarantees
 * a smooth change in the cost function which is nice/needed for gradient-based
 * optimization.
 * 
 * 
 * Using the notation \alpha_k^i to denote the probability in the i-th voxel to 
 * belong to the k-th class (which is a function of the mesh deformation), the
 * cost function is given by
 * 
 *   C = \sum_i \sum_k \alpha_k^i [ \log \sigma_k^2 + ( y_i - \mu_k )^2 / \sigma_k^2 ]
 *       /
 *       \sum_i \sum_k \alpha_k^i
 * 
 * where y_i denotes the intensity in the y-th voxel, and the means and variances are
 * the ones optimizing C, i.e.,
 *
 *   \mu_k = \sum_i \alpha_k^i y_i / \sum_i \alpha_k^i
 *
 *   \sigma_k^2 = \sum_i \alpha_k^i ( y_i - \mu_k )^2 / \sum_i \alpha_k^i
 *  .
 * Note that \sum_k \alpha_k^i can be smaller than one because of the ignore-the-
 * background-class PVI trick mentioned earlier.
 * 
 * Using the notation
 *   
 *   N_k = \sum_i \alpha_k^i  
 * 
 *   L_k = \sum_i \alpha_k^i y_i  

 *   Q_k = \sum_i \alpha_k^i y_i^2
 * 
 * and
 *
 *   N = \sum_k N_k 
 *
 * for the total (non-integer because some voxels don't have full unity weight) "number" 
 * of voxels, and 
 *
 *   E_k = \log [ \sigma_k^2 + 1 ]
 *
 * as the differential entropy of the Gaussian of the k-th class, the cost function can be
 * re-written as
 *
 *   C = \sum_k N_k/N E_k
 * .
 * Noting that 
 *
 *   \mu_k = L_k / N_k
 * 
 * and 
 *  
 *   \sigma_k^2 = ( Q_k - N_k \mu_k^2 ) / N_k = ( Q_k N_k - L_k^2 ) / N_k^2
 * 
 * this can be written as a pure function of the quantities N_k, L_k and Q_k:
 *
 *   C = \sum_k N_k/N [ \log( Q_k N_k - L_k^2 ) - 2 \log N_k + 1 ]
 *
 * and therefore the gradient (with respect to the mesh node positions) gradC
 * can be written as a linear combination of gradL_k, gradN_k, and gradQ_k:
 *
 *   gradC = \sum_k gradN_k/N E_k 
 *           +
 *           -1/N^2 ( \sum_k gradN_k ) \sum_k N_k E_k
 *           +
 *           \sum_k N_k/N [ 1/( Q_k N_k - L_k )^2 ( gradQ_k N_k + Q_k gradN_k - 2 L_k gradL_k ) - 2/N_k gradN_k ]
 *
 * ; grouping terms of gradL_k, gradN_k, and gradQ_k we finally get
 * 
 *   gradC = 1/N [ sum_k gradN_k W_k^N + gradL_k W_k^L + gradQ_k W_k^Q ]
 *
 * with weights
 *
 *   W_k^N = E_k - C + Q_k/\sigma_k^2/N_k - 2
 *
 *   W_k^L = -2L_k / \sigma_k^2/N_k
 *
 *   W_k^Q = 1/sigma_k^2
 *
 * .
 * A straightforward implementation of this class would therefore have a small helper class to compute
 * these weights and cost C first (using a simple rasterization pass of the mesh), after which the main
 * class then does another rasterization pass to compute and combine the relevant gradient components
 * in each tetrahedron. However, this requires two rasterization steps; instead my implementation here
 * only uses one and gradually builds up gradL_k, gradN_k, and gradQ_k simultaneously with their weights
 * by aggregating each voxels' contribution; in the very end (once the final weights are known) gradL_k,
 * gradN_k, and gradQ_k are then linearly combined.
 *
 * Although this seemed like a good idea at the time, the fact that for each class three gradients need
 * to be stored for each of many threads separately, makes the implementation a big mess and also quite 
 * slow (I think); it also breaks the logical flow of other cost functions (implemented in the superclass)
 * so also this adds another level of complexity (and causes even more messy code).
 *
 */
class ConditionalGaussianEntropyCostAndGradientCalculator: public AtlasMeshPositionCostAndGradientCalculator
{
public :
  
  /** Standard class typedefs */
  typedef ConditionalGaussianEntropyCostAndGradientCalculator  Self;
  typedef AtlasMeshPositionCostAndGradientCalculator Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( ConditionalGaussianEntropyCostAndGradientCalculator, AtlasMeshPositionCostAndGradientCalculator );

  /** Some typedefs */
  typedef itk::Image< float, 3 >  ImageType;

  /** */  
  void SetImage( const ImageType* image );

  /** */  
  void Rasterize( const AtlasMesh* mesh );
  
  
protected:
  ConditionalGaussianEntropyCostAndGradientCalculator();
  virtual ~ConditionalGaussianEntropyCostAndGradientCalculator();

  //
  bool RasterizeTetrahedron( const AtlasMesh* mesh, 
                             AtlasMesh::CellIdentifier tetrahedronId,
                             int threadNumber );
  
private:
  ConditionalGaussianEntropyCostAndGradientCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  //
  ImageType::ConstPointer  m_Image;
  
#if 0  
  //
  std::vector< double >  m_Ns;
  std::vector< double >  m_Ls;
  std::vector< double >  m_Qs;
  
  //
  std::vector< AtlasPositionGradientContainerType::Pointer >  m_NGradients;
  std::vector< AtlasPositionGradientContainerType::Pointer >  m_LGradients;
  std::vector< AtlasPositionGradientContainerType::Pointer >  m_QGradients;
#else
  //
  std::vector< std::vector< ThreadAccumDataType > >  m_ThreadSpecificNs;
  std::vector< std::vector< ThreadAccumDataType > >  m_ThreadSpecificLs;
  std::vector< std::vector< ThreadAccumDataType > >  m_ThreadSpecificQs;
  
  //
  std::vector< std::vector< AtlasPositionGradientThreadAccumContainerType::Pointer > >  m_ThreadSpecificNGradients;
  std::vector< std::vector< AtlasPositionGradientThreadAccumContainerType::Pointer > >  m_ThreadSpecificLGradients;
  std::vector< std::vector< AtlasPositionGradientThreadAccumContainerType::Pointer > >  m_ThreadSpecificQGradients;
  
  //
  std::vector< ThreadAccumDataType >  m_ThreadSpecificPriorCosts;
  std::vector< AtlasPositionGradientThreadAccumContainerType::Pointer >  m_ThreadSpecificPriorGradients;

#endif  
  
  
};


} // end namespace kvl

#endif
