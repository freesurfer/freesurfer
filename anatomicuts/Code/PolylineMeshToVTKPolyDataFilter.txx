#include "PolylineMeshToVTKPolyDataFilter.h"

#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkFieldData.h>
#include <vtkSmartPointer.h>


	template <class TMesh>
		PolylineMeshToVTKPolyDataFilter<TMesh>
		::PolylineMeshToVTKPolyDataFilter()
		{
			m_Output = vtkPolyData::New();
			//this->SetNumberOfInputs (1);
			unsigned char color[3] = {219,112,147};
			this->m_color =  color;
		}


	template <class TMesh>
		PolylineMeshToVTKPolyDataFilter<TMesh>
		::~PolylineMeshToVTKPolyDataFilter()
		{
			m_Output->Delete();
		}


	template <class TMesh>
		void
		PolylineMeshToVTKPolyDataFilter<TMesh>
		::SetInput (MeshType *mesh)
		{
			this->SetNthInput (0, mesh);
			this->Modified();
		}


	template <class TMesh>
		void
		PolylineMeshToVTKPolyDataFilter<TMesh>
		::GenerateData()
		{
			if (TMesh::PointDimension!=3)
				itkExceptionMacro (<<"Only meshes of point dimension 3 are supported");

			m_Output->Initialize();
			m_Output->Allocate();

			vtkUnsignedCharArray* allColors = vtkUnsignedCharArray::New();
			allColors->SetNumberOfComponents (3);
			const MeshType *input = static_cast<MeshType*>(this->GetInput(0));

			vtkSmartPointer<vtkPoints> points = vtkPoints::New();

			/*    for (unsigned int i=0; i<input->GetNumberOfPoints(); i++)
			      {
			      PointType pt;
			      pt.Fill (0.0);
			      input->GetPoint (i, &pt);
			      points->InsertNextPoint ( pt[0], pt[1], pt[2] );
			      }
			      */
			typedef typename MeshType::CellsContainer::ConstIterator ConstCellIterator;
			ConstCellIterator itCell = input->GetCells()->Begin();
			while( itCell!=input->GetCells()->End() )
			{
				CellType *cell = itCell.Value();
				vtkIdType *ids = new vtkIdType [cell->GetNumberOfPoints()];
				//      std::cout << " num " << cell->GetNumberOfPoints() << std::endl;
				int index = 0;
				for (typename CellType::PointIdIterator it=cell->PointIdsBegin(); it!=cell->PointIdsEnd(); it++)
				{
					//      std::cout << " id " << *it;
					PointType pt;
					pt.Fill (0.0);
					input->GetPoint (*it, &pt);
					points->InsertPoint (*it, pt[0], pt[1], pt[2] );
					ids[index] = *it;
					index++;
				}
				//      std::cout << std::endl;

				m_Output->InsertNextCell (VTK_POLY_LINE, cell->GetNumberOfPoints(), ids);
				allColors->InsertNextTuple3 ( this->m_color[0], this->m_color[1], this->m_color[2] );
				delete [] ids;

				++itCell;
			}

			vtkSmartPointer<vtkIntArray> intArrayCellData  = vtkIntArray::New();
			typedef typename MeshType::CellDataContainer CellDataContainer;
			if(  input->GetCellData() != 0 )
			{
				//std::cout << " o " << input->GetCellData()->Begin() << std::endl;

				typename CellDataContainer::ConstIterator cellData = input->GetCellData()->Begin();
				int i=0;
				for(;cellData!= input->GetCellData()->End(); ++cellData)
				{
					//      std::cout <<  " cell data " << cellData.Value() << std::endl;
					//  intArrayCellData->InsertNextTuple1(cellData.Value());
					intArrayCellData->InsertValue(i,cellData.Value());
					i++;
				}

				vtkSmartPointer<vtkFieldData> fieldData = vtkFieldData::New();
				fieldData->AddArray(intArrayCellData);
				m_Output->SetFieldData(fieldData);

				//vtkSmartPointer<vtkIntArray> hola = (vtkIntArray*)m_Output->GetFieldData()->GetArray(0);
				//  std::cout << "  get value " <<hola->GetValue(0) << std::endl;
				//int v;
				//hola->GetTupleValue(0,&v) ;
				//std::cout << " get v " << v << std::endl;
			}
			m_Output->SetPoints ( points );
			m_Output->GetCellData()->SetScalars ( allColors );
			points->Delete();

		}



