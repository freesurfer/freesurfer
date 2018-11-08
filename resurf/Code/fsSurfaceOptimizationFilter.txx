#ifndef __fsSurfaceOptimizationFilter_txx
#define __fsSurfaceOptimizationFilter_txx

template< typename TInputSurface , typename TOutputSurface>
void fs::SurfaceOptimizationFilter< TInputSurface, TOutputSurface>::GenerateData()
{
	//this->GetOutput()->SetPoints(this->GetInput()->GetPoints());

	const typename InSurfaceType::PointsContainer* points = this->GetInput()->GetPoints();
	typename InSurfaceType::PointsContainerConstIterator it = points->Begin();
	typename InSurfaceType::PointsContainerConstIterator itEnd = points->End();
	for(unsigned int i=0;it!=itEnd;it++,i++)
	{
		typename InSurfaceType::PointType point = it.Value();
		this->GetOutput()->SetPoint( i, point);
	
	}
	for(int w=0;w<2000;w++)
	{
		points = this->GetOutput()->GetPoints();
		it = points->Begin();
		itEnd = points->End();
		for(unsigned int i=0;it!=itEnd;it++,i++)
		{
			typename InSurfaceType::PointType point = it.Value();
			std::vector<typename InSurfaceType::PointType> adj = this->GetInput()->GetAdjacentPoints(i);	

			typename InSurfaceType::PointType grad = point.EuclideanDistanceTo(adj[0]) - point.EuclideanDistanceTo(adj[1]);
			//std::cout << point.EuclideanDistanceTo(adj[0]) << point.EuclideanDistanceTo(adj[1]) <<std::endl;
			//typename InSurfaceType::PointType temp =  point.GetVectorFromOrigin()*2;
			//temp = temp - adj[0];
			//temp= temp -adj[1];
			//std::cout << adj[0] << adj[1]<< std::endl;
			for ( int j=0;j<3;j++)
			{
				grad[j] = point[j]+ 0.0001*  grad[j] * (2*point[j] - adj[0][j] -adj[1][j]);
			}
			this->GetOutput()->SetPoint( i, grad);

		}
	}
	const typename InSurfaceType::CellsContainer* cells = this->GetInput()->GetCells( );
	typename InSurfaceType::CellsContainerConstIterator itCells = cells->Begin();
	typename InSurfaceType::CellsContainerConstIterator itCellsEnd = cells->End();
	for(unsigned int i=0 ;itCells != itCellsEnd;itCells++,i++ )
	{
		typename InSurfaceType::CellType::CellAutoPointer cellCopy;
		itCells.Value()->MakeCopy( cellCopy );
		this->GetOutput()->SetCell( i, cellCopy );
	}

} 
#endif
