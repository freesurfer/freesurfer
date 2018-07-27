#ifndef _VTKPolyDataToPolylineMeshFilter_txx_
#define _VTKPolyDataToPolylineMeshFilter_txx_

#include "VTKPolyDataToPolylineMeshFilter.h"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkFieldData.h>


  template <class TImage>
  VTKPolyDataToPolylineMeshFilter<TImage>
  ::VTKPolyDataToPolylineMeshFilter()
  {
    this->ProcessObject::SetNumberOfRequiredOutputs(1);
    m_VTKPolyData = 0;
  typename TImage::Pointer outputMesh = TImage::New();
    this->ProcessObject::SetNthOutput(0, outputMesh.GetPointer());

}



  template<class TImage>
  void
  VTKPolyDataToPolylineMeshFilter<TImage>
  ::GenerateData()
  {

	GenerateData2();
    typename OutputMeshType::Pointer outputMesh = this->GetOutput();
     vtkFieldData *fieldData = m_VTKPolyData->GetFieldData();
    int k=0;
    if(fieldData != 0)
    {

        for( int i=0;i<fieldData->GetNumberOfArrays();i++)
        {
            vtkIntArray *arrayCellData = (vtkIntArray*) fieldData->GetArray(i);
            for(int j=0; j< arrayCellData->GetNumberOfTuples();j++)
             {
			std::vector<int> hola;
			hola.push_back(arrayCellData->GetValue(j));
//                outputMesh->SetCellData(k, hola);
                  outputMesh->SetCellData(k, arrayCellData->GetValue(j));

                    k++;
            }
        }
    }
    // cell data array "Labels"
    
       if (vtkDataArray *array = m_VTKPolyData->GetCellData()->GetArray ("Labels"))
       {
       typedef typename OutputMeshType::CellDataContainer CellDataContainer;
       outputMesh->SetCellData( CellDataContainer::New() );
       outputMesh->GetCellData()->Reserve ( outputMesh->GetNumberOfCells() );

       for (int i=0; i<array->GetNumberOfTuples(); i++)
       {
     double val = *(array->GetTuple (i));
     	std::vector<int> hola;
			hola.push_back(*(array->GetTuple (i)));
//                   outputMesh->SetCellData(i, hola);
outputMesh->SetCellData (i, val);
       }
       }
     

  }


  template<class TImage>
  void
  VTKPolyDataToPolylineMeshFilter<TImage>
  ::GenerateData2()
  {

    if (TImage::PointDimension!=3)
      itkExceptionMacro (<<"Only meshes of point dimension 3 are supported");

    if (!m_VTKPolyData)
      itkExceptionMacro (<<"VTK polydata is not set");
   
    typename OutputMeshType::Pointer outputMesh = this->GetOutput();
    outputMesh->SetCellsAllocationMethod(
        OutputMeshType::CellsAllocatedDynamicallyCellByCell );

    outputMesh->GetPoints()->Reserve ( m_VTKPolyData->GetNumberOfPoints() );


    PointIdentifier pointId = 0;
    CellIdentifier  cellId  = 0;

    vtkCellArray *lines = m_VTKPolyData->GetLines();
    lines->InitTraversal();
    vtkIdType pointCount, *pointBuf;
    while ( lines->GetNextCell(pointCount, pointBuf) )
    {

      CellAutoPointer line;
      line.TakeOwnership ( new PolylineCellType );

      for (vtkIdType k=0; k<pointCount; k++)
      {
        double *pt = m_VTKPolyData->GetPoint ( pointBuf[k] );

        PointType point;
        for (unsigned int i=0; i<TImage::PointDimension; i++)
        {
          point [i] = pt[i];
        }

        outputMesh->SetPoint (pointId, point);
        line->SetPointId ( (PointIdentifier)k, pointId);

        pointId++;
      }
      outputMesh->SetCell (cellId++, line);

    }

  }


  template<class TImage>
  void
  VTKPolyDataToPolylineMeshFilter<TImage>
  ::PrintSelf( std::ostream& os, Indent indent ) const
  {
    Superclass::PrintSelf(os,indent);
    
    os << indent << "VTK PolyData: " << *m_VTKPolyData << std::endl;
  }
  


#endif
