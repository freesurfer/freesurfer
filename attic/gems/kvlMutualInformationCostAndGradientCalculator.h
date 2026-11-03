#ifndef __kvlMutualInformationCostAndGradientCalculator_h
#define __kvlMutualInformationCostAndGradientCalculator_h

#include "kvlAtlasMeshPositionCostAndGradientCalculator.h"
#include "kvlHistogrammer.h"


namespace kvl
{


/**
 *
 * Implementation of the D'Agostino MICCAI 2004 flavor of probabilistic-atlas-based Mutual Information, 
 * with some useful tricks/hacks happily borrowed from John Ashburner's spm_maff8.m implementation in
 * SPM12, including:
 * 
 *    1. Iterative estimation of the class-conditional distributions with EM. Standard "partial volume
 *       interpolation" (PVI) based Mutual Information implementations (cf. Maes TMI 1997) effectively only
 *       do a single iteration of this
 * 
 *   2. (Very) approximate gradients that completely ignore the effect of changing the registration position
 *       on the class-conditional distributions; of voxels entering/leaving the overlapping regions in
 *       which the histogram is computed; and even of non-smooth entering/leaving of such voxels (which is
 *       what the PVI thing was originally all about, if I recall correctly) 
 * 
 * The implementation could be sped up by moving the EM stuff into the kvlHistogrammer class, so that 
 * the priors need to be rasterized only once (using existing classes); however I'm really too lazy for 
 * that now.
 * 
 */
class MutualInformationCostAndGradientCalculator: public AtlasMeshPositionCostAndGradientCalculator
{
public :
  
  /** Standard class typedefs */
  typedef MutualInformationCostAndGradientCalculator  Self;
  typedef AtlasMeshPositionCostAndGradientCalculator Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( MutualInformationCostAndGradientCalculator, AtlasMeshPositionCostAndGradientCalculator );

  /** Some typedefs */
  typedef itk::Image< float, 3 >  ImageType;

  /** */  
  void SetImage( const ImageType* image );

  /** */  
  void Rasterize( const AtlasMesh* mesh );
  
  
protected:
  MutualInformationCostAndGradientCalculator();
  virtual ~MutualInformationCostAndGradientCalculator();
  
  void AddDataContributionOfTetrahedron( const AtlasMesh::PointType& p0,
                                         const AtlasMesh::PointType& p1,
                                         const AtlasMesh::PointType& p2,
                                         const AtlasMesh::PointType& p3,
                                         const AtlasAlphasType&  alphasInVertex0,
                                         const AtlasAlphasType&  alphasInVertex1,
                                         const AtlasAlphasType&  alphasInVertex2,
                                         const AtlasAlphasType&  alphasInVertex3,
                                         ThreadAccumDataType&  priorPlusDataCost,
                                         AtlasPositionGradientThreadAccumType&  gradientInVertex0,
                                         AtlasPositionGradientThreadAccumType&  gradientInVertex1,
                                         AtlasPositionGradientThreadAccumType&  gradientInVertex2,
                                         AtlasPositionGradientThreadAccumType&  gradientInVertex3 );
  
private:
  MutualInformationCostAndGradientCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  //
  Histogrammer::Pointer  m_Histogrammer;
  double  m_NumberOfVoxels;
  
};


} // end namespace kvl

#endif
