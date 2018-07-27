#ifndef __VTKPolyDataToPolylineMeshFilter_h
#define __VTKPolyDataToPolylineMeshFilter_h

#include "itkMesh.h"
#include "itkMeshSource.h"
#include "itkPolylineCell.h"

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

  template <class TOutputMesh>
    class VTKPolyDataToPolylineMeshFilter : public MeshSource<TOutputMesh>
  {
  public:
    typedef VTKPolyDataToPolylineMeshFilter Self;
    typedef MeshSource<TOutputMesh>         Superclass;
    typedef SmartPointer<Self>              Pointer;
    typedef SmartPointer<const Self>        ConstPointer;

    itkNewMacro  (Self);
    itkTypeMacro (VTKPolyDataToPolylineMeshFilter, MeshSource);

    typedef TOutputMesh OutputMeshType;
    typedef typename OutputMeshType::MeshTraits MeshTraits;
    typedef typename OutputMeshType::PointType  PointType;
    typedef typename OutputMeshType::PixelType  PixelType;

    /** Some convenient typedefs. */
    typedef typename OutputMeshType::Pointer         OutputMeshPointer;
    typedef typename OutputMeshType::CellTraits      CellTraits;
    typedef typename OutputMeshType::CellIdentifier  CellIdentifier;
    typedef typename OutputMeshType::CellType        CellType;
    typedef typename OutputMeshType::CellAutoPointer CellAutoPointer;
    typedef typename OutputMeshType::PointIdentifier PointIdentifier;
    typedef typename CellTraits::PointIdIterator     PointIdIterator;
    
    typedef typename OutputMeshType::PointsContainerPointer
      PointsContainerPointer;
    
    typedef typename OutputMeshType::PointsContainer
      PointsContainer;

    typedef itk::PolylineCell<CellType>                      PolylineCellType;
    typedef typename  PolylineCellType::SelfAutoPointer SelfAutoPointer;

    void SetVTKPolyData (vtkSmartPointer<vtkPolyData> polydata)
    { m_VTKPolyData = polydata; }
    vtkSmartPointer<vtkPolyData> GetVTKPolyData (void) const
    { return m_VTKPolyData; }
	void GenerateData2();
  protected:
    VTKPolyDataToPolylineMeshFilter();
    ~VTKPolyDataToPolylineMeshFilter() {}
    void PrintSelf(std::ostream& os, Indent indent) const;

    /** Reads the file */
    void GenerateData();

  private:
    VTKPolyDataToPolylineMeshFilter (const Self&);
    void operator=(const Self&);

    vtkSmartPointer<vtkPolyData> m_VTKPolyData;
    
  };

#include "VTKPolyDataToPolylineMeshFilter.txx"  
#endif
