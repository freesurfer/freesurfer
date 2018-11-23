#include "kvlAtlasMeshJacobianDeterminantDrawer.h"

#include "kvlTetrahedronInteriorIterator.h"


namespace kvl
{

//
//
//
AtlasMeshJacobianDeterminantDrawer
::AtlasMeshJacobianDeterminantDrawer()
{
  m_Image = 0; 
}  



//
//
//
AtlasMeshJacobianDeterminantDrawer
::~AtlasMeshJacobianDeterminantDrawer()
{
  
}  

 
//
//
//
bool
AtlasMeshJacobianDeterminantDrawer
::RasterizeTetrahedron( const AtlasMesh* mesh, 
                        AtlasMesh::CellIdentifier tetrahedronId,
                        int threadNumber )
{
  //
  // Retrieve everything we need to know 
  //
  AtlasMesh::CellAutoPointer  cell;
  mesh->GetCell( tetrahedronId, cell );

  AtlasMesh::CellType::PointIdIterator  pit = cell->PointIdsBegin();
  const AtlasMesh::PointIdentifier  id0 = *pit;
  ++pit;
  const AtlasMesh::PointIdentifier  id1 = *pit;
  ++pit;
  const AtlasMesh::PointIdentifier  id2 = *pit;
  ++pit;
  const AtlasMesh::PointIdentifier  id3 = *pit;
  
  AtlasMesh::PointType p0;
  AtlasMesh::PointType p1;
  AtlasMesh::PointType p2;
  AtlasMesh::PointType p3;
  mesh->GetPoint( id0, &p0 );
  mesh->GetPoint( id1, &p1 );
  mesh->GetPoint( id2, &p2 );
  mesh->GetPoint( id3, &p3 );

  const ReferenceTetrahedronInfo&  info = mesh->GetCellData()->ElementAt( tetrahedronId );

  
  //
  // Compute the Jacobian determinant of this tetrahedron
  //
  // Z is inv( [ p0 p1 p2 p3; 1 1 1 1 ] ) of the tetrahedron in reference position
  const double  z11 = info.m_Z11;
  const double  z21 = info.m_Z21;
  const double  z31 = info.m_Z31;
  const double  z41 = info.m_Z41;

  const double  z12 = info.m_Z12;
  const double  z22 = info.m_Z22;
  const double  z32 = info.m_Z32;
  const double  z42 = info.m_Z42;

  const double  z13 = info.m_Z13;
  const double  z23 = info.m_Z23;
  const double  z33 = info.m_Z33;
  const double  z43 = info.m_Z43;

  //  Y = [ p0 p1 p2 p3; 1 1 1 1 ]
  const double y11 = p0[ 0 ];
  const double y21 = p0[ 1 ];
  const double y31 = p0[ 2 ];

  const double y12 = p1[ 0 ];
  const double y22 = p1[ 1 ];
  const double y32 = p1[ 2 ];

  const double y13 = p2[ 0 ];
  const double y23 = p2[ 1 ];
  const double y33 = p2[ 2 ];
  
  const double y14 = p3[ 0 ];
  const double y24 = p3[ 1 ];
  const double y34 = p3[ 2 ];

  // M = Y * Z
  const double  m11 = z11*y11 + z21*y12 + z31*y13 + z41*y14;
  const double  m21 = z11*y21 + z21*y22 + z31*y23 + z41*y24;
  const double  m31 = z11*y31 + z21*y32 + z31*y33 + z41*y34;
  const double  m12 = z12*y11 + z22*y12 + z32*y13 + z42*y14;
  const double  m22 = z12*y21 + z22*y22 + z32*y23 + z42*y24;
  const double  m32 = z12*y31 + z22*y32 + z32*y33 + z42*y34;
  const double  m13 = z13*y11 + z23*y12 + z33*y13 + z43*y14;
  const double  m23 = z13*y21 + z23*y22 + z33*y23 + z43*y24;
  const double  m33 = z13*y31 + z23*y32 + z33*y33 + z43*y34;

  const double  detJ = m11 * ( m22*m33 - m32*m23 ) - m12 * ( m21*m33 - m31*m23 ) + m13 * ( m21*m32 - m31*m22 );

  
  
  // Loop over all voxels within the tetrahedron and do The Right Thing  
  TetrahedronInteriorIterator< ImageType::PixelType >  it( m_Image, p0, p1, p2, p3 );
  for ( ; !it.IsAtEnd(); ++it )
    {
    it.Value() = static_cast< ImageType::PixelType >( detJ );
    }
    
  return true;
}


  
} // End namespace kvl
