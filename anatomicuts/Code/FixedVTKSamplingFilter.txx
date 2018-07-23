#ifndef _itk_FixedVTKSamplingFilter_txx_
#define _itk_FixedVTKSamplingFilter_txx_

#include "interpolation.h"
#include "itkFixedVTKSamplingFilter.h"
#include <vtkCellArray.h>
#include <vtkCellData.h>


namespace itk
{
	template <class TInputMesh, class TOutputMesh>
		FixedVTKSamplingFilter<TInputMesh, TOutputMesh>
		::FixedVTKSamplingFilter()
		{
			this->SetNumberOfRequiredInputs (1);
			this->m_extendPercentage=0.0;	
		}


		template <class TInputMesh, class TOutputMesh>
		void
		FixedVTKSamplingFilter<TInputMesh, TOutputMesh>
		::GenerateData()
		{
			typedef typename TInputMesh::CellsContainer CellsContainer;
			//std::cout << this->GetInput()->GetNumberOfCells() << std::endl;
			this->GetInput()->GetCells()->Begin();
			typename CellsContainer::ConstIterator cellIt = this->GetInput()->GetCells()->Begin();
			//initializing output mesh
			typename OutputMeshType::Pointer outputMesh = this->GetOutput();
			outputMesh->SetCellsAllocationMethod(
					OutputMeshType::CellsAllocatedDynamicallyCellByCell );
			if(this->GetSampling()>0)
			{
				outputMesh->GetPoints()->Reserve ( this->GetInput()->GetNumberOfCells()*this->GetSampling());
			}
			//      std::cout << " number of cells " << this->GetInput()->GetNumberOfCells() << std::endl;
//			std::cout << " sampling " << this->GetSampling() << std::endl;
			OutputPointIdentifier index = 0;
			OutputCellIdentifier indexCell = 0;
			for(;cellIt!= this->GetInput()->GetCells()->End(); ++cellIt)
			{
				int j=0;
				double x,y,z;
				double vecPoints[3*cellIt.Value()->GetNumberOfPoints()];
				typename TInputMesh::CellTraits::PointIdIterator  pointIdIt  = cellIt.Value()->PointIdsBegin();
				alglib::real_2d_array pts;
				typename TInputMesh::PointType ptPrev;  
				this->GetInput()->GetPoint(*pointIdIt, &ptPrev);
				double dist=0;
				for(;pointIdIt != cellIt.Value()->PointIdsEnd();pointIdIt++)
				{
					typename TInputMesh::PointType pt;  
					this->GetInput()->GetPoint(*pointIdIt, &pt);
					for(int i=0;i<3;i++)
					{
						vecPoints[j] =pt[i];
						j++;
					}
					dist += ptPrev.SquaredEuclideanDistanceTo(pt);
					ptPrev = pt;
				}
				pts.setcontent(cellIt.Value()->GetNumberOfPoints(),3, vecPoints);

				alglib::pspline3interpolant s;
				alglib::pspline3build(pts,cellIt.Value()->GetNumberOfPoints(),1,0,s);

				float ii=0;
				int k=0;
				OutputCellAutoPointer line;
				line.TakeOwnership ( new PolylineCellType );
				int sampling = ceil(dist);
				if(this->GetSampling()>1)
				{
					sampling = this->GetSampling();
				}
				else
				{
					sampling = dist / this->GetSampling()*-1;
				}
				//	std::cout <<  "sampling fixed vtk filter" <<  sampling << std::endl;	
				ii=-this->m_extendPercentage;
//				while(k<sampling)
				while(k<sampling)
				{
					alglib::pspline3calc(s,ii,x,y,z);
					typename TInputMesh::PointType pt;
					pt[0]=x;
					pt[1]=y;
					pt[2]=z;
					 //     std::cout << "pt" << pt << std::endl;
					ii+=((1.0 + 2* this->m_extendPercentage)/(sampling-1.));
					outputMesh->SetPoint (index, pt);
					line->SetPointId ( (OutputPointIdentifier)k, index);

					k++;
					index++;

				}
				outputMesh->SetCell (indexCell, line);
				indexCell+=1;
			}    

			outputMesh->SetCellData(const_cast<typename TInputMesh::CellDataContainer*>(this->GetInput()->GetCellData()));

		}

	template <class TInputMesh, class TOutputMesh>
		void
		FixedVTKSamplingFilter<TInputMesh, TOutputMesh>
		::SampleToFiberLenght()
		{
			typedef typename TInputMesh::CellsContainer FixedCellsContainer;
			typename FixedCellsContainer::ConstIterator cellItFxd = this->GetInput()->GetCells()->Begin();
			typename TInputMesh::PointsContainer::Iterator it=this->GetInput()->GetPoints()->Begin(); 
			typename TInputMesh::PointType ptPrevious = it.Value();  

			double fiberLenght = 0;
			it++;
			for(;it != this->GetInput()->GetPoints()->End();it++)
			{
				double dist=0;
				for(int i=0;i<3;i++)
				{
					dist += pow(it.Value()[i] - ptPrevious[i],2);
				}
				dist = sqrt(dist);
				fiberLenght +=dist;
				ptPrevious = it.Value();
			}
			this->SetSampling( ceil(fiberLenght));

		}

}
#endif
