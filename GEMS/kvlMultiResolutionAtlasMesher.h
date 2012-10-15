/**
 * @file  kvlMultiResolutionAtlasMesher.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:40 $
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
#ifndef __kvlMultiResolutionAtlasMesher_h
#define __kvlMultiResolutionAtlasMesher_h

#include "kvlAtlasParameterEstimator.h"


namespace kvl
{


class MultiResolutionAtlasMesher: public itk::Object
{
public :

  /** Standard class typedefs */
  typedef MultiResolutionAtlasMesher  Self;
  typedef itk::Object  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( MultiResolutionAtlasMesher, itk::Object );

  // Some typedefs
  typedef itk::Image< unsigned char, 3 >  LabelImageType;

  // Set label images.
  void SetLabelImages( const std::vector< LabelImageType::ConstPointer >& labelImages );

  // Get label images
  const std::vector< LabelImageType::ConstPointer >&  GetLabelImages() const
  {
    return m_LabelImages;
  }

  // Get label image
  const LabelImageType*  GetLabelImage( unsigned int labelImageNumber ) const;

  // Set/Get initial size
  const unsigned int* GetInitialSize() const
  {
    return m_InitialSize;
  }

  //
  void  SetNumberOfUpsamplingSteps( unsigned int  numberOfUpsamplingSteps )
  {
    m_NumberOfUpsamplingSteps = numberOfUpsamplingSteps;
  }
  unsigned int  GetNumberOfUpsamplingSteps() const
  {
    return m_NumberOfUpsamplingSteps;
  }

  //
  void  SetTryToBeSparse( bool tryToBeSparse )
  {
    m_TryToBeSparse = tryToBeSparse;
  }
  bool  GetTryToBeSparse() const
  {
    return m_TryToBeSparse;
  }

  // Set/Get initial stiffness
  float  GetInitialStiffness( unsigned int upsampleStepNumber ) const
  {
    return m_InitialStiffnesses[ upsampleStepNumber ];
  }

  //
  void  SetUp( const unsigned int* size, const float* initialStiffnesses );

  //
  unsigned int GetNumberOfClasses() const
  {
    return m_NumberOfClasses;
  }

  //
  unsigned int GetNumberOfMeshes() const
  {
    return m_NumberOfMeshes;
  }


  //
  const AtlasMeshCollection*  GetCurrentMeshCollection() const
  {
    return m_Current;
  }

  //
  const AtlasParameterEstimator* GetEstimator() const
  {
    return m_Estimator;
  }

  //
  void SetPositionEstimationIterationEventResolution( unsigned int resolution )
  {
    m_Estimator->SetPositionEstimationIterationEventResolution( resolution );
  }
  unsigned int GetPositionEstimationIterationEventResolution() const
  {
    return m_Estimator->GetPositionEstimationIterationEventResolution();
  }

  //
  void SetPositionGradientDescentStepSize( float stepSize )
  {
    m_Estimator->SetPositionGradientDescentStepSize( stepSize );
  }
  float  GetPositionGradientDescentStepSize() const
  {
    return m_Estimator->GetPositionGradientDescentStepSize();
  }

  //
  void SetAlphaEstimationStopCriterion( float stopCriterion )
  {
    m_Estimator->SetAlphaEstimationStopCriterion( stopCriterion );
  }
  float GetAlphaEstimationStopCriterion() const
  {
    return m_Estimator->GetAlphaEstimationStopCriterion();
  }

  //
  void SetPositionEstimationStopCriterion( float stopCriterion )
  {
    m_Estimator->SetPositionEstimationStopCriterion( stopCriterion );
  }
  float GetPositionEstimationStopCriterion() const
  {
    return m_Estimator->GetPositionEstimationStopCriterion();
  }

  //
  unsigned long AddEstimatorObserver( const itk::EventObject& event, itk::Command* command )
  {
    return m_Estimator->AddObserver( event, command );
  }

  //
  void  Go();

protected :
  // Constructor
  MultiResolutionAtlasMesher();

  // Destructor
  virtual ~MultiResolutionAtlasMesher();

  // Print
  void PrintSelf( std::ostream& os, itk::Indent indent ) const;

  //
  virtual void SetUp();

  //
  typedef itk::AutomaticTopologyMeshSource< kvl::AtlasMesh >  MeshSourceType;
  static void AddHexahedron( MeshSourceType* meshSource,
                             AtlasMesh::CellsContainer* hexahedra,
                             const AtlasMesh::PointType&  p0,
                             const AtlasMesh::PointType&  p1,
                             const AtlasMesh::PointType&  p2,
                             const AtlasMesh::PointType&  p3,
                             const AtlasMesh::PointType&  p4,
                             const AtlasMesh::PointType&  p5,
                             const AtlasMesh::PointType&  p6,
                             const AtlasMesh::PointType&  p7 )
  {
    AtlasMesh::PointIdentifier  p0IdDummy;
    AtlasMesh::PointIdentifier  p1IdDummy;
    AtlasMesh::PointIdentifier  p2IdDummy;
    AtlasMesh::PointIdentifier  p3IdDummy;
    AtlasMesh::PointIdentifier  p4IdDummy;
    AtlasMesh::PointIdentifier  p5IdDummy;
    AtlasMesh::PointIdentifier  p6IdDummy;
    AtlasMesh::PointIdentifier  p7IdDummy;

    AddHexahedron( meshSource, hexahedra,
                   p0, p1, p2, p3, p4, p5, p6, p7,
                   p0IdDummy, p1IdDummy, p2IdDummy, p3IdDummy, p4IdDummy,
                   p5IdDummy, p6IdDummy, p7IdDummy );

  }

  //
  static void AddHexahedron( MeshSourceType* meshSource,
                             AtlasMesh::CellsContainer* hexahedra,
                             const AtlasMesh::PointType&  p0,
                             const AtlasMesh::PointType&  p1,
                             const AtlasMesh::PointType&  p2,
                             const AtlasMesh::PointType&  p3,
                             const AtlasMesh::PointType&  p4,
                             const AtlasMesh::PointType&  p5,
                             const AtlasMesh::PointType&  p6,
                             const AtlasMesh::PointType&  p7,
                             AtlasMesh::PointIdentifier&  p0Id,
                             AtlasMesh::PointIdentifier&  p1Id,
                             AtlasMesh::PointIdentifier&  p2Id,
                             AtlasMesh::PointIdentifier&  p3Id,
                             AtlasMesh::PointIdentifier&  p4Id,
                             AtlasMesh::PointIdentifier&  p5Id,
                             AtlasMesh::PointIdentifier&  p6Id,
                             AtlasMesh::PointIdentifier&  p7Id );




  template<class TCoordRep>
  static void GetUpsampledHexahedronPoints( const AtlasMesh::PointType&  p0,
      const AtlasMesh::PointType&  p1,
      const AtlasMesh::PointType&  p2,
      const AtlasMesh::PointType&  p3,
      const AtlasMesh::PointType&  p4,
      const AtlasMesh::PointType&  p5,
      const AtlasMesh::PointType&  p6,
      const AtlasMesh::PointType&  p7,
      AtlasMesh::PointType&  p01,
      AtlasMesh::PointType&  p02,
      AtlasMesh::PointType&  p13,
      AtlasMesh::PointType&  p23,
      AtlasMesh::PointType&  p15,
      AtlasMesh::PointType&  p37,
      AtlasMesh::PointType&  p26,
      AtlasMesh::PointType&  p04,
      AtlasMesh::PointType&  p57,
      AtlasMesh::PointType&  p67,
      AtlasMesh::PointType&  p46,
      AtlasMesh::PointType&  p45,
      AtlasMesh::PointType&  p0123,
      AtlasMesh::PointType&  p1357,
      AtlasMesh::PointType&  p2367,
      AtlasMesh::PointType&  p0246,
      AtlasMesh::PointType&  p0145,
      AtlasMesh::PointType&  p4567,
      AtlasMesh::PointType&  pMiddle );


private :
  //
  MultiResolutionAtlasMesher(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  //
  AtlasMesh::CellsContainer::Pointer  GetCells( const AtlasMesh::PointsContainer* position ) const;

  //
  AtlasMeshCollection::Pointer  GetMeshCollection( AtlasMesh::PointsContainer* referencePosition,
      std::vector< AtlasMesh::PointsContainer::Pointer >& positions,
      float stiffness ) const;

  //
  void Upsample();

  //
  std::vector< LabelImageType::ConstPointer >  m_LabelImages;
  unsigned int  m_InitialSize[ 3 ];
  unsigned int  m_NumberOfUpsamplingSteps;
  bool  m_TryToBeSparse;
  float  m_InitialStiffnesses[ 5 ];
  AtlasParameterEstimator::Pointer  m_Estimator;

  unsigned int  m_NumberOfClasses;
  unsigned int  m_NumberOfMeshes;
  unsigned int  m_DomainSize[ 3 ];

  AtlasMeshCollection::Pointer  m_Current;

  AtlasMesh::CellsContainer::Pointer  m_Hexahedra;


};



} // end namespace kvl


#endif
