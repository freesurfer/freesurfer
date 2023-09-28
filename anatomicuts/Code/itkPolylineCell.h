/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPolylineCell.h,v $
  Language:  C++
  Date:      $Date: 2010-09-03 13:51:02 -0400 (Fri, 03 Sep 2010) $
  Version:   $Revision: 84 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkPolylineCell_h
#define __itkPolylineCell_h

#include "itkCellInterface.h"
#include "itkVertexCell.h"
#include "itkVectorContainer.h"

namespace itk
{

/** \class PolylineCell
 * PolylineCell represents a line segment for a Mesh.
 *
 * Template parameters for PolylineCell:
 *
 * TPixelType =
 *     The type associated with a point, cell, or boundary for use in storing
 *     its data.
 *
 * TCellTraits =
 *     Type information of mesh containing cell.
 *
 * \ingroup MeshObjects
 */

template < typename TCellInterface >
class ITK_EXPORT PolylineCell: public TCellInterface
{
public:
  /** Standard class typedefs. */
  itkCellCommonTypedefs(PolylineCell);
  itkCellInheritedTypedefs(TCellInterface);
  
  /** Standard part of every itk Object. */
  itkTypeMacro(PolylineCell, CellInterface);

  /** The type of boundary for this lines's vertices. */
  typedef VertexCell< TCellInterface >         VertexType;
  typedef typename VertexType::SelfAutoPointer VertexAutoPointer;

  typedef VectorContainer<unsigned int, PointIdentifier>   PointIdentifierVector;
  typedef typename PointIdentifierVector::Pointer          PointIdentifierVectorPointer;
    
  /** Polyline-specific topology numbers. */
  enum { CellDimension = 1 };

  
  /** Implement the standard CellInterface. */
#if ITK_VERSION_MAJOR >= 5 && ITK_VERSION_MINOR >= 4
  virtual itk::CommonEnums::CellGeometry GetType() const
#else  
  virtual CellGeometry GetType(void) const
#endif    
  {return Superclass::POLYGON_CELL;}
  virtual void MakeCopy( CellAutoPointer & ) const;
  virtual unsigned int GetDimension(void) const;
  virtual unsigned int GetNumberOfPoints(void) const;
  virtual CellFeatureCount GetNumberOfBoundaryFeatures(int dimension) const;
  virtual bool GetBoundaryFeature(int dimension, CellFeatureIdentifier,CellAutoPointer &);
  virtual void SetPointIds(PointIdConstIterator first);
  virtual void SetPointIds(PointIdConstIterator first,
                           PointIdConstIterator last);
  virtual void SetPointId(int localId, PointIdentifier);
  virtual PointIdIterator      PointIdsBegin(void);
  virtual PointIdConstIterator PointIdsBegin(void) const;
  virtual PointIdIterator      PointIdsEnd(void);
  virtual PointIdConstIterator PointIdsEnd(void) const; 
  
  /** Line-specific interface. */
  virtual CellFeatureCount GetNumberOfVertices(void) const;
  virtual bool GetVertex(CellFeatureIdentifier, VertexAutoPointer &);
  
  /** Visitor interface */
  itkCellVisitMacro(Superclass::POLYGON_CELL);

  PolylineCell()
  {
    m_PointIds = PointIdentifierVector::New();
  }
  ~PolylineCell() {}

protected:
  /** Store number of points needed for a line segment. */
  PointIdentifierVectorPointer m_PointIds;

private:
  PolylineCell(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};


} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPolylineCell.txx"
#endif

#endif
