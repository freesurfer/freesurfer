#ifndef _itk_FilterFibers_txx_
#define _itk_FilterFibers_txx_
#include <random>
#include <chrono>
#include "FilterFibers.h"
#include "itkPolylineCell.h"

template <class TInputMesh, class TOutputMesh>
void
	FilterFibers<TInputMesh, TOutputMesh>
::GenerateData()
{
	const InputMeshType *input = this->GetInput();

	typename InputCellsContainer::ConstIterator inputCellIt = input->GetCells()->Begin();
	int cellId =0;
	typename OutputMeshType::Pointer outputMesh = this->GetOutput();
	outputMesh->SetCellsAllocationMethod( OutputMeshType::CellsAllocatedDynamicallyCellByCell );
	typedef std::mt19937 generator;
	//	typedef std::tr1::mt19937 eng;
	//			std::default_random_engine generator;
	//			std::uniform_real_distribution<float> uniform(0, this->GetInput()->GetNumberOfCells());
	std::uniform_int<double> uniform(0,this->GetInput()->GetNumberOfCells());
	//	 std::tr1::uniform_int<int> unif(1, 52);


	std::mt19937 eng(std::chrono::high_resolution_clock::now()
			.time_since_epoch().count());
	std::uniform_int<int> unif(0,this->GetInput()->GetNumberOfCells());
	std::cout << std::endl << this->GetInput()->GetNumberOfCells() << " " << unif(eng) << '\n';
	int index=0;
	while(cellId < this->GetPercentage()*this->GetInput()->GetNumberOfCells())
	{

		int number = unif(eng) % this->GetInput()->GetNumberOfCells();
		//	int number = unif();
		//		std::cout << " number " << number << std::endl;
		typename OutputMeshType::CellAutoPointer cell;
		typename OutputMeshType::CellAutoPointer line;
		typedef itk::PolylineCell<typename OutputMeshType::CellType>    PolylineCellType;
		line.TakeOwnership ( new PolylineCellType);

		this->GetInput()->GetCell(number, cell);		
		int k=0;
		for (typename OutputCellType::PointIdIterator it = cell->PointIdsBegin(); it!=cell->PointIdsEnd(); it++)
		{
			InputPointType pt;
			pt.Fill (0.0);
			input->GetPoint (*it, &pt);

			outputMesh->SetPoint (index,pt);
			line->SetPointId ( (typename OutputMeshType::PointIdentifier)k, index);
			k++;
			index++;
		}    

		outputMesh->SetCell(cellId,line);
		cellId++;
	}

}


#endif
