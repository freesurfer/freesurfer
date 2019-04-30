#include "itkNumericTraitsVariableLengthVectorPixel.h"
#include "vnl/vnl_math.h"
#include <cstring>
#include <stdlib.h>
#include <string.h>


template< typename TValueType ,class TMesh >
void LabelPerPointVariableLengthVector<TValueType,TMesh>::SetCell( MeshPointerType mesh,  int cellId) 
{
	int pointBegin =0;
	int pointEnd= mesh->GetNumberOfPoints();
	if( cellId != -1)
	{
		CellAutoPointerType cellAutoPointer;
		mesh->GetCell(cellId, cellAutoPointer);	
		////typename MeshType::CellTraits::PointIdIterator  pointIdIt  = cellAutoPointer->PointIdsBegin();
		pointBegin=cellId*cellAutoPointer->GetNumberOfPoints();
		pointEnd= (cellId+1)*cellAutoPointer->GetNumberOfPoints();
	}
	int i=0;
	typename MeshType::PointType ptPrevious=0;  
	mesh->GetPoint(pointBegin, &ptPrevious);
	this->m_length =0.0;
	this->m_numberOfPoints = pointEnd-pointBegin;
	//std::cout << this->m_numberOfPoints<< std::endl;
		for(int pointIdIt =pointBegin; pointIdIt < pointEnd;pointIdIt++)
	{

		typename MeshType::PointType pt=0;  
		mesh->GetPoint(pointIdIt, &pt);
		this->m_length += ptPrevious.EuclideanDistanceTo(pt);
		ptPrevious = pt;
		(*this)[i]= pt[0];
		(*this)[i+1]= pt[1];
		(*this)[i+2]= pt[2];

		i+=3;
		CellType cell;
		if(mesh->GetPointData()->IndexExists(pointIdIt))
		{
			CellType cell = mesh->GetPointData()->GetElement(pointIdIt);
			this->m_labels.push_back(cell);

			if( this->m_labelsPerDirection.size() < cell.size())
			{
				for(int q=0;q<cell.size();q++)
					this->m_labelsPerDirection.push_back(LabelsMapType()); 
			}
			for(int q=0;q<cell.size();q++)
			{
				if (this->m_labelsPerDirection[q].count( cell[q] )>0)
				{
					this->m_labelsPerDirection[q][cell[q]]= this->m_labelsPerDirection[q][cell[q]]+1.0; //000.0/this->m_numberOfPoints;
					//std::cout << this->m_labelsPerDirection[q][cell[q]] << std::endl;
				}
				else
				{
					this->m_labelsPerDirection[q][cell[q]]= 1.0; //000.0/this->m_numberOfPoints;
					//std::cout <<1.0/this->m_numberOfPoints << " " << this->m_labelsPerDirection[q][cell[q]] << std::endl;
				}
			}
		}
	}
	this->m_cellId = cellId;

}
template< typename TValueType ,class TMesh >
void LabelPerPointVariableLengthVector<TValueType,TMesh>
::Print() const
{
	std::cout << " holaaaa " ;
	for(int i=0;i<this->m_labels.size();i++)
	{
		for(int j=0;j<this->m_labels[i].size();j++)
			std::cout <<  this->m_labels[i][j] << " ";
	}
	std::cout << std::endl;

}




