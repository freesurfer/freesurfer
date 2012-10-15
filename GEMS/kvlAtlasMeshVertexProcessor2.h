/**
 * @file  kvlAtlasMeshVertexProcessor2.h
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
#ifndef __kvlAtlasMeshVertexProcessor_h
#define __kvlAtlasMeshVertexProcessor_h

#include "kvlAtlasMeshCollection.h"
#include "itkImage.h"



namespace kvl
{


class AtlasMeshVertexProcessor: public itk::Object
{
public :

  /** Standard class typedefs */
  typedef AtlasMeshVertexProcessor  Self;
  typedef itk::Object  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshVertexProcessor, itk::Object );

  // Some typedefs
  typedef itk::Image< unsigned char, 3 >  LabelImageType;

  //
  void SetMeshCollection( AtlasMeshCollection* meshCollection )
  {
    m_MeshCollection = meshCollection;
  }

  AtlasMeshCollection* GetMeshCollection()
  {
    return m_MeshCollection;
  }

  const AtlasMeshCollection* GetMeshCollection() const
  {
    return m_MeshCollection;
  }

  // Set label images.
  void SetLabelImages( const std::vector< LabelImageType::ConstPointer >& labelImages );

  // Get label images
  const std::vector< LabelImageType::ConstPointer >&  GetLabelImages() const
  {
    return m_LabelImages;
  }

  // Set beta
  void  SetBeta( float beta )
  {
    m_Beta = beta;
  }

  // Get beta
  float  GetBeta() const
  {
    return m_Beta;
  }

  //
  bool CalculateXstar( int meshNumber, AtlasMesh::PointIdentifier pointId,
                       float& xstar, float& ystar, float& zstar, bool verbose=false ) const;

  //
  float  CalculateCost( int meshNumber, AtlasMesh::PointIdentifier pointId, float x, float y, float z ) const;

  //
  AtlasPositionGradientType  CalculateGradient( int meshNumber, AtlasMesh::PointIdentifier pointId, float x, float y, float z ) const;

  //
  float  CalculateCostAndGradient( int meshNumber, AtlasMesh::PointIdentifier pointId, float x, float y, float z,
                                   AtlasPositionGradientType& gradient ) const;
  //
  Curvature  CalculateCurvature( int meshNumber, AtlasMesh::PointIdentifier pointId, float x, float y, float z, bool verbose=false ) const;


protected :
  // Constructor
  AtlasMeshVertexProcessor();

  // Destructor
  virtual ~AtlasMeshVertexProcessor();

  // Print
  void PrintSelf( std::ostream& os, itk::Indent indent ) const;

#if 0
  //
  Curvature  CalculatePriorCurvature( int meshNumber, AtlasMesh::PointIdentifier pointId, float x, float y, float z );
#endif


private :
  AtlasMeshVertexProcessor(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  // Data members
  AtlasMeshCollection::Pointer  m_MeshCollection;
  float  m_Beta;
  std::vector< LabelImageType::ConstPointer >  m_LabelImages;


};



} // end namespace kvl


#endif
