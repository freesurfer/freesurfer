#include "kvlAtlasMeshValueDrawer.h"
#include "kvlTetrahedronInteriorIterator.h"


namespace kvl {


bool AtlasMeshValueDrawer::RasterizeTetrahedron(const AtlasMesh* mesh, AtlasMesh::CellIdentifier tetrahedronId, int threadNumber)
{
  // get tetrahedron cell
  AtlasMesh::CellAutoPointer cell;
  mesh->GetCell(tetrahedronId, cell);

  // get point IDs
  AtlasMesh::CellType::PointIdIterator pit = cell->PointIdsBegin();
  const AtlasMesh::PointIdentifier id0 = *pit;
  ++pit;
  const AtlasMesh::PointIdentifier id1 = *pit;
  ++pit;
  const AtlasMesh::PointIdentifier id2 = *pit;
  ++pit;
  const AtlasMesh::PointIdentifier id3 = *pit;

  // get points
  AtlasMesh::PointType p0;
  AtlasMesh::PointType p1;
  AtlasMesh::PointType p2;
  AtlasMesh::PointType p3;
  mesh->GetPoint(id0, &p0);
  mesh->GetPoint(id1, &p1);
  mesh->GetPoint(id2, &p2);
  mesh->GetPoint(id3, &p3);
  
  const double * valuesInVertex0 = &m_Values[id0 * m_NumFrames];
  const double * valuesInVertex1 = &m_Values[id1 * m_NumFrames];
  const double * valuesInVertex2 = &m_Values[id2 * m_NumFrames];
  const double * valuesInVertex3 = &m_Values[id3 * m_NumFrames];

  // loop over all voxels within the tetrahedron and interpolate

  // ----
  TetrahedronInteriorIterator<ImageType::PixelType> it(m_Image, p0, p1, p2, p3);
  for (int f = 0; f < m_NumFrames; f++) 
    {
    if (valuesInVertex0[ f ] != 0 || valuesInVertex1[ f ] != 0 || valuesInVertex2[ f ] != 0 || valuesInVertex3[ f ] != 0)
      it.AddExtraLoading(valuesInVertex0[f], valuesInVertex1[f], valuesInVertex2[f], valuesInVertex3[f]);
    }
  
  // ----
  for ( ; !it.IsAtEnd(); ++it ) {
    it.Value() = AtlasAlphasType(m_NumFrames);
    int fIdx = 0;
    for (int f = 0; f < m_NumFrames; f++)
      { 
      if (valuesInVertex0[ f ] != 0 || valuesInVertex1[ f ] != 0 || valuesInVertex2[ f ] != 0 || valuesInVertex3[ f ] != 0)
	{
        it.Value()[f] = it.GetExtraLoadingInterpolatedValue(fIdx);
        fIdx++;
	}
      else
        it.Value()[f] = 0;
      }
  }
    
  return true;
}

  
} // End namespace kvl
