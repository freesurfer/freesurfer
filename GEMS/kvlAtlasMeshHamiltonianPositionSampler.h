/**
 * @file  kvlAtlasMeshHamiltonianPositionSampler.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:38 $
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
#ifndef __kvlAtlasMeshHamiltonianPositionSampler_h
#define __kvlAtlasMeshHamiltonianPositionSampler_h

#include "kvlAtlasMeshCollection.h"
#include "itkImage.h"


namespace kvl
{


class AtlasMeshHamiltonianPositionSampler: public itk::Object
{
public :

  /** Standard class typedefs */
  typedef AtlasMeshHamiltonianPositionSampler  Self;
  typedef itk::Object  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshHamiltonianPositionSampler, itk::Object );

  // Some typedefs
  typedef itk::Image< unsigned char, 3 >  ImageType;

  //
  void SetMeshCollection( AtlasMeshCollection* meshCollection )
  {
    m_MeshCollection = meshCollection;
  }

  //
  const AtlasMeshCollection* GetMeshCollection() const
  {
    return m_MeshCollection;
  }

  // Set image
  void SetImage( const ImageType*  image )
  {
    m_Image = image;
  }

  // Get image
  const ImageType*  GetImage() const
  {
    return m_Image;
  }

  //
  void SetMeans( const std::vector< float >&  means )
  {
    m_Means = means;
  }
  const std::vector< float >&  GetMeans() const
  {
    return m_Means;
  }

  //
  void SetVariances( const std::vector< float >&  variances )
  {
    m_Variances = variances;
  }
  const std::vector< float >&  GetVariances() const
  {
    return m_Variances;
  }

  //
  void SetNumberOfWarmupSweeps( unsigned int numberOfWarmupSweeps )
  {
    m_NumberOfWarmupSweeps = numberOfWarmupSweeps;
  }
  unsigned int GetNumberOfWarmupSweeps() const
  {
    return m_NumberOfWarmupSweeps;
  }

  //
  void SetNumberOfBetweenSamplesSweeps( unsigned int numberOfBetweenSamplesSweeps )
  {
    m_NumberOfBetweenSamplesSweeps = numberOfBetweenSamplesSweeps;
  }
  unsigned int GetNumberOfBetweenSamplesSweeps() const
  {
    return m_NumberOfBetweenSamplesSweeps;
  }

  //
  void SetNumberOfSamples( unsigned int numberOfSamples )
  {
    m_NumberOfSamples = numberOfSamples;
  }
  unsigned int GetNumberOfSamples() const
  {
    return m_NumberOfSamples;
  }

  //
  void Reseed( int seed );

  //
  void  SetInitialPositionNumber( int initialPositionNumber )
  {
    m_InitialPositionNumber = initialPositionNumber;
  }

  //
  int  GetInitialPositionNumber() const
  {
    return m_InitialPositionNumber;
  }

  //
  void  SetMaximalTrackingTime( float maximalTrackingTime )
  {
    m_MaximalTrackingTime = maximalTrackingTime;
  }

  //
  float GetMaximalTrackingTime() const
  {
    return m_MaximalTrackingTime;
  }

  //
  void SetTimeStep( float timeStep )
  {
    m_TimeStep = timeStep;
  }

  //
  float  GetTimeStep() const
  {
    return m_TimeStep;
  }

  //
  void  SetMass( float mass )
  {
    m_Mass = mass;
  }

  //
  float  GetMass() const
  {
    return m_Mass;
  }

  // Let the beast go
  AtlasMeshCollection::Pointer  GetSamples();


protected :
  // Constructor
  AtlasMeshHamiltonianPositionSampler();

  // Destructor
  virtual ~AtlasMeshHamiltonianPositionSampler();

  // Print
  void PrintSelf( std::ostream& os, itk::Indent indent ) const;

  //
  AtlasPositionGradientContainerType::Pointer  GetGradient( const AtlasMesh::PointsContainer*  position,
      double& cost ) const;


private :
  AtlasMeshHamiltonianPositionSampler(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  // Data members
  AtlasMeshCollection::Pointer  m_MeshCollection;
  ImageType::ConstPointer  m_Image;
  std::vector< float >  m_Means;
  std::vector< float >  m_Variances;

  unsigned int m_NumberOfWarmupSweeps;
  unsigned int m_NumberOfBetweenSamplesSweeps;
  unsigned int m_NumberOfSamples;
  int  m_InitialPositionNumber;
  float  m_MaximalTrackingTime;
  float  m_TimeStep;
  float  m_Mass;

};



} // end namespace kvl


#endif
