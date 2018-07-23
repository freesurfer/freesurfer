/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPolylineCell.txx,v $
  Language:  C++
  Date:      $Date: 2010-09-03 13:51:02 -0400 (Fri, 03 Sep 2010) $
  Version:   $Revision: 84 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkPolylineCell_txx
#define __itkPolylineCell_txx
#include "itkPolylineCell.h"

namespace itk
{

/**
 * Standard CellInterface:
 */
template <typename TCellInterface>
void
PolylineCell< TCellInterface >
::MakeCopy(CellAutoPointer & cellPointer) const
{
  cellPointer.TakeOwnership( new Self );
  cellPointer->SetPointIds(this->PointIdsBegin(), this->PointIdsEnd());
}


  
/**
 * Standard CellInterface:
 * Get the topological dimension of this cell.
 */
template <typename TCellInterface>
unsigned int
PolylineCell< TCellInterface >
::GetDimension(void) const
{
  return Self::CellDimension;
}

/**
 * Standard CellInterface:
 * Get the number of points required to define the cell.
 */
template <typename TCellInterface>
unsigned int
PolylineCell< TCellInterface >
::GetNumberOfPoints(void) const
{
  return m_PointIds->Size();
}  


/**
 * Standard CellInterface:
 * Get the number of boundary entities of the given dimension.
 */
template <typename TCellInterface>
typename PolylineCell< TCellInterface >::CellFeatureCount
PolylineCell< TCellInterface >
::GetNumberOfBoundaryFeatures(int dimension) const
{
  switch (dimension)
    {
    case 0: return GetNumberOfVertices();
    default: return 0;
    }
}


/**
 * Standard CellInterface:
 * Get the boundary feature of the given dimension specified by the given
 * cell feature Id.
 * The Id can range from 0 to GetNumberOfBoundaryFeatures(dimension)-1.
 */
template <typename TCellInterface>
bool
PolylineCell< TCellInterface >
::GetBoundaryFeature(int dimension, CellFeatureIdentifier featureId, 
                     CellAutoPointer & cellPointer)
{
  switch (dimension)
    {
    case 0: 
      {
      VertexAutoPointer vertexPointer;
      if( this->GetVertex(featureId,vertexPointer) )
        {
        TransferAutoPointer(cellPointer,vertexPointer);
        return true;
        }
      else
        {
        cellPointer.Reset();
        return false;
        }
      break;
      }
    default: 
      {
      cellPointer.Reset();
      return false;
      }
    }
  return false;
}


/**
 * Standard CellInterface:
 * Set the point id list used by the cell.  It is assumed that the given
 * iterator can be incremented and safely de-referenced enough times to 
 * get all the point ids needed by the cell.
 */
template <typename TCellInterface>
void
PolylineCell< TCellInterface >
::SetPointIds(PointIdConstIterator first)
{
  PointIdConstIterator ii(first);
  for(unsigned int i=0; i < m_PointIds->Size(); ++i)
    {
      m_PointIds->CreateIndex (i);
      m_PointIds->SetElement (i, *ii++);
    }
}


/**
 * Standard CellInterface:
 * Set the point id list used by the cell.  It is assumed that the range
 * of iterators [first, last) contains the correct number of points needed to
 * define the cell.  The position *last is NOT referenced, so it can safely
 * be one beyond the end of an array or other container.
 */
template <typename TCellInterface>
void
PolylineCell< TCellInterface >
::SetPointIds(PointIdConstIterator first, PointIdConstIterator last)
{
  int localId=0;
  PointIdConstIterator ii(first);

  m_PointIds->Initialize();
  while(ii != last)
    {
      m_PointIds->CreateIndex (localId);
      m_PointIds->SetElement (localId++, *ii++);
    }
}


/**
 * Standard CellInterface:
 * Set an individual point identifier in the cell.
 */
template <typename TCellInterface>
void
PolylineCell< TCellInterface >
::SetPointId(int localId, PointIdentifier ptId)
{
  m_PointIds->CreateIndex (localId);
  m_PointIds->SetElement (localId, ptId);
}


/**
 * Standard CellInterface:
 * Get a begin iterator to the list of point identifiers used by the cell.
 */
template <typename TCellInterface>
typename PolylineCell< TCellInterface >::PointIdIterator
PolylineCell< TCellInterface >
::PointIdsBegin(void)
{
  return &m_PointIds->ElementAt (0);
}


/**
 * Standard CellInterface:
 * Get a const begin iterator to the list of point identifiers used
 * by the cell.
 */
template <typename TCellInterface>
typename PolylineCell< TCellInterface >::PointIdConstIterator
PolylineCell< TCellInterface >
::PointIdsBegin(void) const
{
  return &m_PointIds->ElementAt (0);
}


/**
 * Standard CellInterface:
 * Get an end iterator to the list of point identifiers used by the cell.
 */
template <typename TCellInterface>
typename PolylineCell< TCellInterface >::PointIdIterator
PolylineCell< TCellInterface >
::PointIdsEnd(void)
{
  return &m_PointIds->ElementAt ( m_PointIds->Size()-1 ) + 1;
}


/**
 * Standard CellInterface:
 * Get a const end iterator to the list of point identifiers used
 * by the cell.
 */
template <typename TCellInterface>
typename PolylineCell< TCellInterface >::PointIdConstIterator
PolylineCell< TCellInterface >
::PointIdsEnd(void) const
{
  return &m_PointIds->ElementAt ( m_PointIds->Size()-1 ) + 1;
}


/**
 * Line-specific:
 * Get the number of vertices for this line.
 */
template <typename TCellInterface>
typename PolylineCell< TCellInterface >::CellFeatureCount
PolylineCell< TCellInterface >
::GetNumberOfVertices(void) const
{
  return this->GetNumberOfPoints();
}


/**
 * Line-specific:
 * Get the vertex specified by the given cell feature Id.
 * The Id can range from 0 to GetNumberOfVertices()-1.
 */
template <typename TCellInterface>
bool
PolylineCell< TCellInterface >
::GetVertex(CellFeatureIdentifier vertexId, VertexAutoPointer & vertexPointer )
{
  VertexType * vert = new VertexType;
  vert->SetPointId(0, m_PointIds->GetElement (vertexId));
  vertexPointer.TakeOwnership( vert );
  return true;  
}

} // end namespace itk

#endif
