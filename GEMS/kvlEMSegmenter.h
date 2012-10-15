/**
 * @file  kvlEMSegmenter.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:39 $
 *    $Revision: 1.3 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */
#ifndef __kvlEMSegmenter_h
#define __kvlEMSegmenter_h

#include "itkImage.h"
#include <vector>
#include "kvlAtlasMesh.h"
#include "kvlImageSmoother.h"


namespace kvl
{



/**
 *
 */
class EMSegmenter: public itk::Object
{
public :

  /** Standard class typedefs */
  typedef EMSegmenter  Self;
  typedef itk::Object  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( EMSegmenter, itk::Object );

  /** Some typedefs */
  typedef itk::Image< unsigned short, 3 >  ImageType;
  typedef itk::Image< float, 3 >  ClassificationImageType;
  typedef itk::Image< float, 3 >  BiasCorrectedImageType;
  typedef itk::Image< float, 3 >  BiasFieldImageType;


  // Set image
  void SetImage( const ImageType*  image );

  // Get image
  const ImageType*  GetImage() const
  {
    return m_Image;
  }

  // Get image
  const BiasCorrectedImageType*  GetBiasCorrectedImage() const
  {
    return m_BiasCorrectedImage;
  }

  // Get estimated bias field image
  const BiasFieldImageType*  GetBiasField() const
  {
    return m_BiasField;
  }


  // Set atlas mesh
  void SetAtlasMesh( const AtlasMesh* atlasMesh,
                     const std::vector< unsigned int >& lookupTable = std::vector< unsigned int >(),
                     const std::vector< unsigned int >& independentParametersLookupTable = std::vector< unsigned int >() );

  // Get atlas mesh
  const AtlasMesh* GetAtlasMesh() const
  {
    return m_AtlasMesh;
  }
  AtlasMesh* GetAtlasMesh()
  {
    return m_AtlasMesh;
  }

  //
  double  Segment();

  //
  void  SetImageAndSegmentItWithCurrentParametersAndWithoutBiasFieldCorrection( const ImageType* image );

  //
  BiasCorrectedImageType::Pointer  BiasCorrect( const ImageType* image, int upsamplingFactor ) const;

  //
  static ImageType::Pointer  ConvertToNativeImageType( const BiasCorrectedImageType* image );

  //
  unsigned int GetNumberOfClasses() const
  {
    return m_NumberOfClasses;
  }

  //
  const ClassificationImageType* GetPrior( unsigned int classNumber ) const
  {
    if ( m_Priors.size() == 0 )
    {
      this->Initialize();
    }
    return m_Priors[ classNumber ];
  }

  //
  const ClassificationImageType* GetPosterior( unsigned int classNumber ) const
  {
    if ( m_Priors.size() == 0 )
    {
      this->Initialize();
    }
    return m_Posteriors[ classNumber ];
  }

  //
  const std::vector< float >&  GetMeans() const
  {
    return m_Means;
  }

  //
  const std::vector< float >&  GetVariances() const
  {
    return m_Variances;
  }

  //
  const std::vector< float >&  GetPriorWeights() const
  {
    return m_PriorWeights;
  }

  //
  int GetIterationNumber() const
  {
    return m_IterationNumber;
  }

  //
  void SetMaximumNumberOfIterations( int maximumNumberOfIterations )
  {
    m_MaximumNumberOfIterations = maximumNumberOfIterations;
  }

  //
  int GetMaximumNumberOfIterations() const
  {
    return m_MaximumNumberOfIterations;
  }

  //
  void  SetStopCriterion( float stopCriterion )
  {
    m_StopCriterion = stopCriterion;
  }

  //
  float  GetStopCriterion() const
  {
    return m_StopCriterion;
  }

  //
  const std::vector< unsigned int >&  GetLookupTable() const
  {
    return m_LookupTable;
  }

  //
  const std::vector< unsigned int >&  GetIndependentParametersLookupTable() const
  {
    return m_IndependentParametersLookupTable;
  }

  //
  void  SetPartialVolumeUpsamplingFactors( const int* partialVolumeUpsamplingFactors )
  {
    m_NumberOfSubvoxels = 1;
    for ( int i = 0; i < 3; i++ )
    {
      m_PartialVolumeUpsamplingFactors[ i ] = partialVolumeUpsamplingFactors[ i ];
      m_NumberOfSubvoxels *= m_PartialVolumeUpsamplingFactors[ i ];

      if ( m_PartialVolumeUpsamplingFactors[ i ] != 1 )
      {
        this->SetMaximumNumberOfIterations( 100 );
      }
    }
  }

  //
  const int*  GetPartialVolumeUpsamplingFactors() const
  {
    return m_PartialVolumeUpsamplingFactors;
  }

  //
  void  SetBiasFieldOrder( int biasFieldOrder )
  {
    m_BiasFieldOrder = biasFieldOrder;
    m_BiasFieldResidueSmoother = 0;
  }

  //
  int  GetBiasFieldOrder() const
  {
    return m_BiasFieldOrder;
  }


  //
  void  SetReestimatedRelativeWeightOfSharedClasses( bool reestimatedRelativeWeightOfSharedClasses )
  {
    m_ReestimatedRelativeWeightOfSharedClasses = reestimatedRelativeWeightOfSharedClasses;
  }

  bool  GetReestimatedRelativeWeightOfSharedClasses() const
  {
    return m_ReestimatedRelativeWeightOfSharedClasses;
  }

  //
  typedef itk::Image< AtlasAlphasType, 3 >  SuperResolutionProbabilityImageType;
  const SuperResolutionProbabilityImageType*  GetSuperResolutionPriors() const
  {
    return m_SuperResolutionPriors;
  }
  SuperResolutionProbabilityImageType::Pointer  GetSuperResolutionLikelihoods() const;
  SuperResolutionProbabilityImageType::Pointer  GetSuperResolutionPosteriors() const;

  //
  ImageType::Pointer  GetSuperResolutionImage() const;

  //
  AtlasMesh::Pointer  GetSuperResolutionMesh( const AtlasMesh* mesh ) const;
  AtlasMesh::Pointer  GetNormalResolutionMesh( const AtlasMesh* superResolutionMesh ) const;


protected:
  EMSegmenter();
  virtual ~EMSegmenter();

  //
  void Initialize() const;

  //
  double UpdatePosteriors();

  //
  void UpdateModelParameters( int biasFieldOrderUsed );

  //
  void UpdateMixtureModelParameters();

  //
  void UpdateBiasFieldModelParameters( int biasFieldOrderUsed );

  //
  static itk::Array< float >  GetMixingProbabilities( const itk::Array< float >& pureProbabilitiesOfClass1,
      const itk::Array< float >& pureProbabilitiesOfClass2 );

  //
  float  DoOneEMIterationWithPartialVoluming();

  //
  static AtlasMesh::Pointer  GetTransformedMesh( const AtlasMesh* mesh, const float* translation, const float* scaling );


private:
  EMSegmenter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  //
  AtlasMesh::Pointer  GetDistributedMesh( const AtlasMesh* atlasMesh, const std::vector< unsigned int >& lookupTable ) const;

  //
  void  SplitSharedClasses();

  // Data members
  ImageType::ConstPointer  m_Image;
  BiasCorrectedImageType::Pointer  m_BiasCorrectedImage;
  BiasFieldImageType::Pointer  m_BiasField;
  AtlasMesh::Pointer  m_AtlasMesh;

  std::vector< unsigned int >  m_LookupTable;
  std::vector< unsigned int >  m_IndependentParametersLookupTable;

  unsigned int  m_NumberOfClasses;
  int  m_BiasFieldOrder;
  bool  m_ReestimatedRelativeWeightOfSharedClasses;

  mutable std::vector< ClassificationImageType::Pointer >  m_Priors;
  mutable std::vector< ClassificationImageType::Pointer >  m_Posteriors;
  std::vector< float >  m_Means;
  std::vector< float >  m_Variances;
  std::vector< float >  m_PriorWeights;

  int m_IterationNumber;
  int m_MaximumNumberOfIterations;
  float  m_StopCriterion;

  // Stuff related to partial volume modeling
  int  m_NumberOfSubvoxels;
  int  m_PartialVolumeUpsamplingFactors[ 3 ];

  mutable SuperResolutionProbabilityImageType::Pointer  m_SuperResolutionPriors;

  // For bias field correction, we have a smoother that fits a polynomial
  // model to a noisy residue image. Initializing the smoother (in particular
  // setting up orthonormal basis functions) takes some time, so we keep it
  // here as a remember that needs to be initialized once
  ImageSmoother::Pointer  m_BiasFieldResidueSmoother;


  typedef itk::Image< std::vector< AtlasAlphasType* >, 3 >  SubvoxelPriorsImageType;
  mutable SubvoxelPriorsImageType::Pointer  m_SubvoxelPriorsImage;  // Low-res image

};


} // end namespace kvl

#endif

