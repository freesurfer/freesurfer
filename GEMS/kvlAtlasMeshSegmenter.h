/**
 * @file  kvlAtlasMeshSegmenter.h
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
#ifndef __kvlAtlasMeshSegmenter_h
#define __kvlAtlasMeshSegmenter_h

#include "kvlEMSegmenter.h"
#include "kvlAtlasMeshToIntensityImagePartialVolumeGradientCalculator.h"
#include "itkAffineTransform.h"



namespace kvl
{

itkEventMacro( PositionUpdatingStartEvent, itk::UserEvent );
itkEventMacro( PositionUpdatingIterationEvent, itk::UserEvent );
itkEventMacro( PositionUpdatingEndEvent, itk::UserEvent );


/**
 *
 */
class AtlasMeshSegmenter: public itk::Object
{
public :

  /** Standard class typedefs */
  typedef AtlasMeshSegmenter  Self;
  typedef itk::Object  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshSegmenter, itk::Object );

  /** Some typedefs */
  typedef itk::Image< unsigned short, 3 >  ImageType;
  typedef itk::Image< float, 3 >  ClassificationImageType;
  typedef itk::AffineTransform< float, 3 >  TransformType;

  //
  void  SetPartialVolumeUpsamplingFactors( const int* partialVolumeUpsamplingFactors );

  //
  const int*  GetPartialVolumeUpsamplingFactors() const
  {
    return m_PartialVolumeUpsamplingFactors;
  }


  // Set image
  void SetImage( const ImageType*  image );

  // Get image
  const ImageType*  GetImage() const
  {
    return m_Image;
  }

  // Set mesh
  void SetMesh( AtlasMesh* mesh,
                const std::vector< unsigned int >&  lookupTable = std::vector< unsigned int >(),
                const std::vector< unsigned int >&  independentParametersLookupTable = std::vector< unsigned int >() );

  // Get atlas mesh
  const AtlasMesh* GetMesh() const
  {
    return m_Mesh;
  }

  //
  unsigned long AddEMSegmenterObserver( const itk::EventObject& event, itk::Command* command )
  {
    return m_EMSegmenter->AddObserver( event, command );
  }

  //
  void Segment( bool useAffine=false );

  //
  unsigned int GetNumberOfClasses() const
  {
    return m_NumberOfClasses;
  }

  //
  const ClassificationImageType* GetPosterior( unsigned int classNumber ) const
  {
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
  const EMSegmenter*  GetEMSegmenter() const
  {
    return m_EMSegmenter;
  }

  //
  unsigned int  GetPositionUpdatingIterationNumber() const
  {
    return m_PositionUpdatingIterationNumber;
  }

  //
  unsigned int  GetPositionUpdatingMaximumNumberOfIterations() const
  {
    return m_PositionUpdatingMaximumNumberOfIterations;
  }

  void  SetPositionUpdatingMaximumNumberOfIterations( unsigned int positionUpdatingMaximumNumberOfIterations )
  {
    m_PositionUpdatingMaximumNumberOfIterations = positionUpdatingMaximumNumberOfIterations;
  }

  //
  unsigned int  GetPositionUpdatingIterationEventResolution() const
  {
    return m_PositionUpdatingIterationEventResolution;
  }

  //
  void SetPositionUpdatingIterationEventResolution( unsigned int  PositionUpdatingIterationEventResolution )
  {
    m_PositionUpdatingIterationEventResolution = PositionUpdatingIterationEventResolution;
  }

  //
  unsigned int  GetIterationNumber() const
  {
    return m_IterationNumber;
  }

  //
  void  SetMaximumNumberOfIterations( unsigned int maximumNumberOfIterations )
  {
    m_MaximumNumberOfIterations = maximumNumberOfIterations;
  }

  unsigned int  GetMaximumNumberOfIterations() const
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
  void  SetMaximalDeformationStopCriterion( float maximalDeformationStopCriterion )
  {
    m_MaximalDeformationStopCriterion = maximalDeformationStopCriterion;
  }

  //
  float  GetMaximalDeformationStopCriterion() const
  {
    return m_MaximalDeformationStopCriterion;
  }

  //
  typedef itk::Image< float, 3 >  SummaryImageType;
  SummaryImageType::Pointer  GetSummaryImage( bool usePrior = false,
      bool useCrisp = false,
      const std::vector< float >&  means = std::vector< float >() ) const;

  //
  void  SetCoregisterToPosteriorProbabilities( bool coregisterToPosteriorProbabilities )
  {
    m_CoregisterToPosteriorProbabilities = coregisterToPosteriorProbabilities;
  }

  bool  GetCoregisterToPosteriorProbabilities() const
  {
    return m_CoregisterToPosteriorProbabilities;
  }

  //
  void  SetBiasFieldOrder( int biasFieldOrder )
  {
    m_EMSegmenter->SetBiasFieldOrder( biasFieldOrder );
  }

  int  GetBiasFieldOrder() const
  {
    return m_EMSegmenter->GetBiasFieldOrder();
  }


  //
  void  SetMeshToImageTranform( const TransformType* meshToImageTransform )
  {
    m_MeshToImageTransform = meshToImageTransform;
  }

  const TransformType* GetMeshToImageTranform() const
  {
    return m_MeshToImageTransform;
  }

  bool  GetDontDeform() const
  {
    return m_DontDeform;
  }

  void  SetDontDeform( bool dontDeform )
  {
    m_DontDeform = dontDeform;
  }

protected:
  AtlasMeshSegmenter();
  virtual ~AtlasMeshSegmenter();

  //
  bool UpdatePosition( bool useAffine );

  //
  float  UpdateModelParameters();

  //
  float  CalculateCurrentWarpCost() const;

  //
  void  MakeAffine( AtlasPositionGradientContainerType* gradient ) const;


  //
  typedef itk::Image< AtlasAlphasType, 3 >  ProbabilityImageType;
  ProbabilityImageType::Pointer  GetPosteriorProbabilityImage() const;


private:
  AtlasMeshSegmenter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  // Data members
  int  m_PartialVolumeUpsamplingFactors[ 3 ];
  bool  m_UsePartialVoluming;

  ImageType::ConstPointer  m_Image;
  typedef itk::Image< unsigned char, 3 >  InternalImageType;
  InternalImageType::ConstPointer  m_InternalImage;
  AtlasMesh::Pointer  m_Mesh;

  unsigned int  m_NumberOfClasses;

  std::vector< ClassificationImageType::ConstPointer >  m_Posteriors;
  std::vector< float >  m_Means;
  std::vector< float >  m_Variances;

  EMSegmenter::Pointer  m_EMSegmenter;

  unsigned int  m_PositionUpdatingIterationNumber;
  unsigned int  m_PositionUpdatingMaximumNumberOfIterations;
  unsigned int  m_PositionUpdatingIterationEventResolution;

  unsigned int  m_IterationNumber;
  unsigned int  m_MaximumNumberOfIterations;
  float  m_StopCriterion;
  float  m_MaximalDeformationStopCriterion;

  bool  m_CoregisterToPosteriorProbabilities;

  TransformType::ConstPointer  m_MeshToImageTransform;

  mutable AtlasMeshToIntensityImagePartialVolumeGradientCalculator::Pointer  m_PartialVolumeGradientCalculator;
  mutable ProbabilityImageType::Pointer  m_PosteriorProbabilityImage;

  bool  m_DontDeform;

};


} // end namespace kvl

#endif

