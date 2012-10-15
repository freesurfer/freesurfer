/**
 * @file  kvlAtlasMeshSegmentationDriver.h
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
#ifndef __kvlAtlasMeshSegmentationDriver_h
#define __kvlAtlasMeshSegmentationDriver_h

#include "kvlAtlasMeshSegmenter.h"
#include "kvlCompressionLookupTable.h"
#include "kvlCroppedImageReader.h"
#include <string>



namespace kvl
{

itkEventMacro( MultiResolutionStartEvent, itk::UserEvent );
itkEventMacro( MultiResolutionIterationEvent, itk::UserEvent );
itkEventMacro( MultiResolutionEndEvent, itk::UserEvent );


/**
 *
 */
class AtlasMeshSegmentationDriver: public itk::Object
{
public :

  /** Standard class typedefs */
  typedef AtlasMeshSegmentationDriver  Self;
  typedef itk::Object  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshSegmentationDriver, itk::Object );

  //
  void  SetUp( const std::string& setUpFileName )
  {
    this->ParseSetUpFile( setUpFileName );
    this->Initialize();
  }

  //
  void Segment( bool  useAffine=false );

  //
  const AtlasMeshSegmenter*  GetAtlasMeshSegmenter() const
  {
    return m_Segmenter;
  }

  //
  const std::vector< std::string >&  GetLabels() const
  {
    return m_Labels;
  }


  //
  unsigned long AddSegmenterObserver( const itk::EventObject& event, itk::Command* command )
  {
    return m_Segmenter->AddObserver( event, command );
  }

  //
  unsigned long AddEMSegmenterObserver( const itk::EventObject& event, itk::Command* command )
  {
    return m_Segmenter->AddEMSegmenterObserver( event, command );
  }


  //
  typedef itk::Image< float, 3 >  SummaryImageType;
  SummaryImageType::Pointer  GetSummaryImage( bool usePrior = false,
      bool useCrisp = false,
      bool joinSplitClasses = true,
      const std::vector< float >&  means = std::vector< float >() ) const;

  //
  void HandleSegmenterEvent( itk::Object* object, const itk::EventObject & event );

protected:
  AtlasMeshSegmentationDriver();
  virtual ~AtlasMeshSegmentationDriver();

  typedef itk::Image< short, 3 >  FreeSurferLabelImageType;
  FreeSurferLabelImageType::Pointer  MapCrispSummaryBackToFreeSurfer( const SummaryImageType*  crispSummary ) const;

  //
  virtual void ParseSetUpFile( const std::string& setUpFileName );

  //
  virtual void Initialize();

  //
  typedef std::map< unsigned char, unsigned char >  MoreToLessLabelsLookupTableType;
  AtlasMesh::ConstPointer  GetDistributedMesh( const AtlasMesh*  originalMesh,
      const MoreToLessLabelsLookupTableType& lookupTable ) const;

  //
  AtlasMesh::ConstPointer  GetReducedMesh( const AtlasMesh*  originalMesh,
      const MoreToLessLabelsLookupTableType& lookupTable ) const;


private:
  AtlasMeshSegmentationDriver(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  // Data members
  AtlasMeshSegmenter::Pointer  m_Segmenter;

  // Stuff to be read from a set up input file
  std::string  m_LogDirectory;
  std::string  m_ImageFileName;
  std::string  m_MeshCollectionFileName;
  std::string  m_BoundingFileName;
  std::string  m_CompressionLookupTableFileName;
  bool  m_CoregisterToPosteriorProbabilities;
  int  m_PartialVolumeUpsamplingFactors[ 3 ];
  int  m_DownSamplingFactor;
  float  m_BackgroundPriorThresholdSmoothingSigma;
  bool  m_LogTransformData;
  int  m_BiasFieldOrder;
  int  m_StartingMeshNumber;
  float  m_K;
  typedef std::map< unsigned int, int > SplitClassesContainerType;
  SplitClassesContainerType  m_SplitClasses;
  typedef std::vector< std::vector< unsigned int > >  SameGaussianParametersContainterType;
  SameGaussianParametersContainterType  m_SameGaussianParameters;
  typedef std::vector< float >  SmoothingSigmasContainterType;
  SmoothingSigmasContainterType  m_MeshSmoothingSigmas;
  SmoothingSigmasContainterType  m_ImageSmoothingSigmas;


  unsigned int  m_MultiResolutionLevel;
  unsigned int  m_NumberOfMultiResolutionLevels;

  typedef AtlasMeshSegmenter::ImageType  ImageType;
  ImageType::Pointer  m_NonDownsampledImage;

  CompressionLookupTable::Pointer  m_Compressor;
  std::vector< std::string >  m_Labels;

  //
  MoreToLessLabelsLookupTableType  m_InternalToReducedLookupTable;
  MoreToLessLabelsLookupTableType  m_InternalToExternalLookupTable;

  ImageType::RegionType  m_OriginalImageCropRegion;
  ImageType::RegionType  m_OriginalImageOriginalRegion;

  //
  typedef CroppedImageReader::TransformType  TransformType;
  TransformType::Pointer  m_WorldToImageTransform;

  //
  static void GetSpacingOriginAndDirection( const TransformType* transform,
      ImageType::SpacingType& spacing,
      ImageType::PointType& origin,
      ImageType::DirectionType& direction );

  // Keep a collection around that can be updated and written out as the internal segmenter
  // loops over its iterations
  AtlasMeshCollection::Pointer  m_TrackingCollection;

};


} // end namespace kvl

#endif

