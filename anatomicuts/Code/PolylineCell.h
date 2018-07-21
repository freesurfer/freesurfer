#include "itkCellInterface.h"
#include "itkVertexCell.h"
#include "itkVectorContainer.h"
using namespace itk;

template < typename TCellInterface >
class  PolylineCell: public TCellInterface
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
  virtual CellGeometry GetType(void) const 
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



