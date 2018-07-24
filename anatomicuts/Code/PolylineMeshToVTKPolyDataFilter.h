#ifndef __PolylineMeshToVTKPolyDataFilter_h
#define __PolylineMeshToVTKPolyDataFilter_h

#include "itkProcessObject.h"
#include "itkPolylineCell.h"
#include <vtkSmartPointer.h>

#include "vtkCellData.h"
using namespace itk;
class vtkPolyData;

  template <class TMesh>
    class PolylineMeshToVTKPolyDataFilter : public ProcessObject
  {
  public:
    typedef PolylineMeshToVTKPolyDataFilter Self;
    typedef ProcessObject                   Superclass;
    typedef SmartPointer<Self>       Pointer;
    typedef SmartPointer<const Self> ConstPointer;

    itkNewMacro  (Self);
    itkTypeMacro (PolylineMeshToVTKPolyDataFilter, ProcessObject);

    typedef TMesh    MeshType;
    typedef typename MeshType::MeshTraits MeshTraits;
    typedef typename MeshType::PointType  PointType;
    typedef typename MeshType::PixelType  PixelType;

    /** Some convenient typedefs. */
    typedef typename MeshType::Pointer         MeshPointer;
    typedef typename MeshType::CellTraits      CellTraits;
    typedef typename MeshType::CellIdentifier  CellIdentifier;
    typedef typename MeshType::CellType        CellType;
    typedef typename MeshType::CellAutoPointer CellAutoPointer;
    typedef typename MeshType::PointIdentifier PointIdentifier;
    typedef typename CellTraits::PointIdIterator     PointIdIterator;
    
    typedef typename MeshType::PointsContainerPointer
      PointsContainerPointer;
    
    typedef typename MeshType::PointsContainer
      PointsContainer;

    typedef PolylineCell<CellType>                      PolylineCellType;
    typedef typename  PolylineCellType::SelfAutoPointer SelfAutoPointer;

    void SetInput (MeshType *mesh);
    
    vtkSmartPointer<vtkPolyData> GetOutputPolyData (void)
    { return m_Output; }

    virtual void Update (void)
    { this->GenerateData(); }
	void SetColor ( unsigned char color[])
	{
		this->m_color = color;
	}

  protected:
    PolylineMeshToVTKPolyDataFilter();
    ~PolylineMeshToVTKPolyDataFilter();

    virtual void GenerateData (void);

  private:
    PolylineMeshToVTKPolyDataFilter (const Self&);
    void operator= (const Self&);
	unsigned char* m_color;
    vtkSmartPointer<vtkPolyData> m_Output;
  };

#include "PolylineMeshToVTKPolyDataFilter.txx"
#endif
