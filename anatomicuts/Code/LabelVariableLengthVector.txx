#ifndef __itkLabelVariableLengthVectorCurrents_txx
#define __itkLabelVariableLengthVectorCurrents_txx

#include "itkVariableLengthVectorCurrents.h"
#include "itkNumericTraitsVariableLengthVectorPixel.h"
#include "vnl/vnl_math.h"
#include <cstring>
#include <stdlib.h>
#include <string.h>


	template< typename TValueType ,class TMesh >
void LabelVariableLengthVector<TValueType,TMesh>::SetCell(MeshPointerType mesh, int cellId)
{
	CellAutoPointerType cellAutoPointer;
	mesh->GetCell(cellId, cellAutoPointer);
	
	//CellType labels;

	typedef typename  MeshType::CellPixelType CellType;
	//std::vector<int> labels;
	CellType labels;
	mesh->GetCellData(cellId, &labels);
	this->SetSize(labels.size());
	for(unsigned int i=0;i<labels.size();i++)
	{
		(*this)[i]=labels[i];
	}	
}


#endif
