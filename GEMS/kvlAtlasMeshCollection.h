/**
 * @file  kvlAtlasMeshCollection.h
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
#ifndef __kvlAtlasMeshCollection_h
#define __kvlAtlasMeshCollection_h

#include "kvlAtlasMesh.h"
#include "itkAffineTransform.h"
#include "itkAutomaticTopologyMeshSource.h"



namespace kvl
{


class AtlasMeshCollection: public itk::Object
{
public :

  /** Standard class typedefs */
  typedef AtlasMeshCollection  Self;
  typedef itk::Object  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshCollection, itk::Object );

  // Some typedefs
  typedef AtlasMesh::PointsContainer  PointsContainerType;
  typedef AtlasMesh::PointDataContainer  PointDataContainerType;
  typedef AtlasMesh::CellsContainer  CellsContainerType;
  typedef AtlasMesh::CellDataContainer  CellDataContainerType;

  //
  void SetPointParameters( PointDataContainerType* pointParameters )
  {
    m_PointParameters = pointParameters;
    m_Meshes.clear();
  }
  const PointDataContainerType*  GetPointParameters() const
  {
    return m_PointParameters;
  }
  PointDataContainerType*  GetPointParameters()
  {
    return m_PointParameters;
  }

  //
  void SetCells( CellsContainerType* cells )
  {
    m_Cells = cells;
    m_Meshes.clear();
    m_CellLinks = 0;
  }
  const CellsContainerType*  GetCells() const
  {
    return m_Cells;
  }
  CellsContainerType*  GetCells()
  {
    return m_Cells;
  }


  //
  void SetReferencePosition( PointsContainerType* referencePosition )
  {
    m_ReferencePosition = referencePosition;
    m_ReferenceTetrahedronInfos = 0;
    m_Meshes.clear();
  }
  const PointsContainerType*  GetReferencePosition() const
  {
    return m_ReferencePosition;
  }
  PointsContainerType*  GetReferencePosition()
  {
    return m_ReferencePosition;
  }

  //
  void SetK( float K )
  {
    m_K = K;
    m_ReferenceTetrahedronInfos = 0;
    m_Meshes.clear();
  }
  float GetK() const
  {
    return m_K;
  }

  //
  const CellDataContainerType*  GetReferenceTetrahedronInfos() const;

  //
  void SetPositions( const std::vector< PointsContainerType::Pointer >& positions )
  {
    m_Positions = positions;
    m_Meshes.clear();
  }
  const std::vector< PointsContainerType::Pointer >&  GetPositions() const
  {
    return m_Positions;
  }
  std::vector< PointsContainerType::Pointer >&  GetPositions()
  {
    return m_Positions;
  }


  //
  unsigned int GetNumberOfMeshes() const
  {
    return static_cast< unsigned int >( m_Positions.size() );
  }

  //
  void GenerateFromSingleMesh( AtlasMesh* mesh, unsigned int numberOfMeshes, float K );

  //
  const AtlasMesh*  GetMesh( unsigned int meshNumber ) const;

  //
  AtlasMesh::ConstPointer  GetReferenceMesh() const;


  // Write out to file
  bool Write( const char* fileName ) const;

  // Read from file
  bool Read( const char* fileName );


  //
  void Construct( const unsigned int*  meshSize, const unsigned int*  domainSize,
                  float K,
                  unsigned int numberOfClasses, unsigned int numberOfMeshes,
                  bool  forceBorderVerticesToBackground = true );


  //
  bool GetCollapsed( AtlasMesh::CellIdentifier edgeId,
                     AtlasMeshCollection::Pointer& collapsed,
                     std::set< AtlasMesh::CellIdentifier >&  disappearingCells,
                     AtlasMesh::CellIdentifier& unifiedVertexId,
                     bool initializeAlphasToFlat = true ) const;


  //
  AtlasMeshCollection::Pointer  GetRegionGrown( AtlasMesh::CellIdentifier seedId,
      unsigned int radius,
      bool makeOuterPointsImmobile = true ) const;

  //
  AtlasMeshCollection::Pointer  GetUpsampled() const;

  //
  AtlasMeshCollection::Pointer  GetEdgeSplitted( AtlasMesh::CellIdentifier  edgeId,
      AtlasMesh::CellIdentifier  newVertexId,
      AtlasMesh::PointIdentifier  newPointId ) const;

  //
  AtlasMeshCollection::Pointer GetEdgeSwapped( AtlasMesh::CellIdentifier edgeId ) const;

  //
  AtlasMeshCollection::Pointer GetEdgeSwapped( AtlasMesh::CellIdentifier edgeId,
      AtlasMesh::CellIdentifier  newVertexId,
      AtlasMesh::PointIdentifier  newPointId,
      AtlasMesh::CellIdentifier&  newEdgeId ) const;


  //
  void  FlattenAlphas();

  //
  typedef  itk::AffineTransform< float, 3 >  TransformType;
  void  Transform( int meshNumber, const TransformType* transform );

  //
  AtlasMesh::CellLinksContainerPointer  GetCellLinks() const;

protected :
  // Constructor
  AtlasMeshCollection();

  // Destructor
  virtual ~AtlasMeshCollection();

  // Print
  void PrintSelf( std::ostream& os, itk::Indent indent ) const;

  //
  AtlasMeshCollection::Pointer  GetEdgeSplitted( AtlasMesh::CellIdentifier  edgeId,
      AtlasMesh::CellIdentifier  newVertexId,
      AtlasMesh::PointIdentifier  newPointId,
      AtlasMesh::CellIdentifier&  transverseEdge0Id,
      AtlasMesh::CellIdentifier&  transverseEdge1Id ) const;

  //
  typedef itk::AutomaticTopologyMeshSource< kvl::AtlasMesh >  MeshSourceType;
  void FillCubeWithTetrahedra( MeshSourceType* meshSource, bool flippedConfiguration,
                               const AtlasMesh::PointType&  p0,
                               const AtlasMesh::PointType&  p1,
                               const AtlasMesh::PointType&  p2,
                               const AtlasMesh::PointType&  p3,
                               const AtlasMesh::PointType&  p4,
                               const AtlasMesh::PointType&  p5,
                               const AtlasMesh::PointType&  p6,
                               const AtlasMesh::PointType&  p7 ) const;

private :
  AtlasMeshCollection(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  // Data members
  PointDataContainerType::Pointer  m_PointParameters;
  std::vector< PointsContainerType::Pointer >  m_Positions;

  CellsContainerType::Pointer  m_Cells;
  PointsContainerType::Pointer m_ReferencePosition;
  float  m_K;
  mutable CellDataContainerType::Pointer  m_ReferenceTetrahedronInfos; // Cached
  mutable std::vector< AtlasMesh::Pointer >  m_Meshes;  // Cached
  mutable AtlasMesh::CellLinksContainerPointer  m_CellLinks;
};



} // end namespace kvl


#endif
