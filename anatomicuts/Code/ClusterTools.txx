#ifndef _ClusterTools_txx_
#define _ClusterTools_txx_

static 	std::string LEFT = "Left";
static	std::string RIGHT = "Right";
static 	std::string LEFT2 = "lh";
static	std::string RIGHT2 = "rh";

COLOR_TABLE *ct=NULL;

template <class TColorMesh, class TImage, class THistogramMesh>
int ClusterTools<TColorMesh, TImage, THistogramMesh>::SymmetricLabelId(int id)
{

	if( ct == NULL)
	{
		FSENV *fsenv = FSENVgetenv();
		char tmpstr[2000];	
		sprintf(tmpstr, "%s/FreeSurferColorLUT.txt", fsenv->FREESURFER_HOME);
		ct = CTABreadASCII(tmpstr);
	}

	std::string str = std::string(ct->entries[id]->name);
	int symId = id;
	if( str.find(LEFT) != std::string::npos)
	{
		char* hola= (char*)str.replace(str.find(LEFT),LEFT.length(), RIGHT).c_str();
		symId = CTABentryNameToIndex(hola, ct); 
	}else if( str.find(LEFT2) != std::string::npos)
	{
		char* hola= (char*)str.replace(str.find(LEFT2),LEFT2.length(), RIGHT2).c_str();
		symId = CTABentryNameToIndex(hola, ct); 

	} else	if( str.find(RIGHT) != std::string::npos)
	{
		char* hola= (char*)str.replace(str.find(RIGHT),RIGHT.length(), LEFT).c_str();
		symId = CTABentryNameToIndex(hola, ct); 
	}else if( str.find(RIGHT2) != std::string::npos)
	{
		char* hola= (char*)str.replace(str.find(RIGHT2),RIGHT2.length(), LEFT2).c_str();
		symId = CTABentryNameToIndex(hola, ct); 

	}
	return  symId;
}
template <class TColorMesh, class TImage, class THistogramMesh>
std::vector<typename THistogramMesh::Pointer> ClusterTools<TColorMesh, TImage, THistogramMesh>::ColorMeshToHistogramMesh(std::vector<typename TColorMesh::Pointer> basicMeshes,typename  TImage::Pointer segmentation, bool removeInterHemispheric)
{
	if( ct == NULL)
	{
		FSENV *fsenv = FSENVgetenv();
		char tmpstr[2000];	
		sprintf(tmpstr, "%s/FreeSurferColorLUT.txt", fsenv->FREESURFER_HOME);
		ct = CTABreadASCII(tmpstr);
	}

	std::vector<typename HistogramMeshType::Pointer> meshes;
	for(unsigned int i=0;i<basicMeshes.size();i++)
	{
		typedef typename ColorMeshType::CellsContainer::ConstIterator CellIterator;
		int globalIndex=0;
		int indexCell =0;
		// typedef HistogramMeshType::PointIdentifier PointIdentifier;
		typedef typename  HistogramMeshType::PointDataContainer PointDataContainerType;
		typename ColorMeshType::Pointer basicMesh = basicMeshes[i];
		typename HistogramMeshType::Pointer mesh = HistogramMeshType::New();
		int out =0;
		for(CellIterator cellIt = basicMesh->GetCells()->Begin(); cellIt!= basicMesh->GetCells()->End(); cellIt++)
		{
			typename PointDataContainerType::Pointer dataContainer = PointDataContainerType::New();
			int withinIndex=0;
			typename HistogramMeshType::CellAutoPointer line;
			typename ColorMeshType::CellTraits::PointIdIterator  pointIdIt  = cellIt.Value()->PointIdsBegin();

			typename ColorMeshType::PointType pt=0; 
			int left=0, right=0;

			if(removeInterHemispheric)
			{ 
				for(pointIdIt  = cellIt.Value()->PointIdsBegin();pointIdIt != cellIt.Value()->PointIdsEnd();pointIdIt++)
				{
					basicMesh->GetPoint(*pointIdIt, &pt);
					typename ImageType::IndexType index;
					if (segmentation->TransformPhysicalPointToIndex(pt,index))
					{
						typename ImageType::PixelType label = segmentation->GetPixel(index);

						std::string str = std::string(ct->entries[(int)label]->name);
						if( str.find(LEFT) != std::string::npos ||  str.find(LEFT2) != std::string::npos)
						{
							left++;						

						} else	if( str.find(RIGHT) != std::string::npos ||  str.find(RIGHT2) != std::string::npos)
						{
							right++;
						}				
					}

				}
			}
			if (right == 0 || left ==0 || !removeInterHemispheric )
			{
				for(pointIdIt  = cellIt.Value()->PointIdsBegin();pointIdIt != cellIt.Value()->PointIdsEnd();pointIdIt++)
				{
					basicMesh->GetPoint(*pointIdIt, &pt);

					line.TakeOwnership ( new itk::PolylineCell<typename HistogramMeshType::CellType> );
					mesh->SetPoint (globalIndex, pt);
					line->SetPointId (withinIndex, globalIndex);

					withinIndex++;
					globalIndex++;
					//std::cout << globalIndex << " "<< std::endl;
				}

				mesh->SetCell (indexCell, line);
				indexCell++;
			}
			else
			{
				out++;
				//	std::cout << "right " << right << " left " << left << " remove "<< removeInterHemispheric << std::endl;
			}
		}
		float val =(float)out / ((float)indexCell+out);
		//std::cout << " val " << val << std::endl;
		if( val> 0.20)
		{
			meshes.push_back(HistogramMeshType::New());
		}
		else
		{
			meshes.push_back(mesh);
		}
	}
	return meshes;
}
/*
template <class TColorMesh, class TImage, class THistogramMesh>
std::vector<MeasurementVectorType> ClusterTools<TColorMesh, TImage, THistogramMesh>::SetDirectionalNeighbors(std::vector<ColorMeshType::Pointer> meshes, std::vector<int> clusterCentroidsIndex, ImageType::Pointer segmentation, std::vector<itk::Vector<float>> direcciones, bool symmetry)
{
	//std::cout << " symmetry " << symmetry << std::endl;
	std::vector<MeasurementVectorType> measurements;
	for(unsigned int i=0;i<meshes.size();i++)
	{	
		int  pointId =0;
		//ColorMeshType::CellTraits::PointIdIterator  pointIdEnd =meshes[i]->GetNumberOfPoints();
		int  pointIdEnd =meshes[i]->GetNumberOfPoints();
		int numPoints = meshes[i]->GetNumberOfPoints();
		if(clusterCentroidsIndex[i]>-1)
		{

			CellAutoPointer cell1;
			meshes[i]->GetCell(clusterCentroidsIndex[i], cell1);
			pointId  =*cell1->PointIdsBegin();
			pointIdEnd=*cell1->PointIdsEnd();
			numPoints= pointIdEnd - pointId;
		}

		for(;pointId != pointIdEnd; pointId++)
		{
			ColorMeshType::PointType pt1=0;
			meshes[i]->GetPoint (pointId, &pt1);
			IndexType index;	
			if (segmentation->TransformPhysicalPointToIndex(pt1,index))
			{
				PointDataType* pointData = new PointDataType();
				PixelType labelOrig = segmentation->GetPixel(index);
				PixelType label = labelOrig;
				if (symmetry)	
				{
					label= SymmetricLabelId(labelOrig);
					//					std::cout << " labeled " << labelOrig  << " mirrowed label " << label << std::endl;
				} 

				pointData->push_back(label);


				for(unsigned int k=0;k<direcciones.size();k++)
				{
					PixelType vecino = labelOrig;
					itk::ContinuousIndex<float,3> continuousIndex = index;
					ColorMeshType::PointType point = pt1;
					//std::cout << direcciones[k] << std::endl;
					IndexType roundedIndex;
					while(vecino == labelOrig)
					{
						for(unsigned int j=0; j<3;j++)
							continuousIndex[j] += direcciones[k][j];
						roundedIndex.CopyWithRound(continuousIndex);
						if(!segmentation->GetLargestPossibleRegion().IsInside(roundedIndex))
							break;
						vecino = segmentation->GetPixel(roundedIndex);
					}
					if(vecino!=0)
					{
						if (symmetry)
						{
							vecino= SymmetricLabelId(vecino);
						}
						pointData->push_back(vecino);
					}
					else
					{
						pointData->push_back(label);
					}
				}
				meshes[i]->SetPointData(pointId, *pointData);
			}
			else
			{

				std::cout << pointId << " " << pointIdEnd << " "<< pt1 <<" " << index<<std::endl;
			}
		}
		MeasurementVectorType mv(numPoints*3);
		mv.SetCell(meshes[i],clusterCentroidsIndex[i]);
		measurements.push_back(mv);

	}
	return measurements;

}
template <class TColorMesh, class TImage, class THistogramMesh>
std::vector<int> ClusterTools<TColorMesh, TImage, THistogramMesh>::GetCentroidIndices(std::vector<ColorMeshType::Pointer> meshes)
{
	std::vector<int> clusterCentroidsIndex( meshes.size(),0);
	for(unsigned int i=0; i<meshes.size();i++)	
	{
		std::vector<ColorMeshType::PointType> avgPoints(meshes[i]->GetCells()->Begin().Value()->GetNumberOfPoints(),0);
		ColorMeshType::CellsContainer::ConstIterator cells = meshes[i]->GetCells()->Begin();
		int cell_i=0;
		for(;cells!= meshes[i]->GetCells()->End(); cells++)
		{
			cell_i++;
			ColorMeshType::CellTraits::PointIdIterator  pointIdIt  =cells.Value()->PointIdsBegin();
			double dist=0.0;
			double dist_inv=0.0;
			for(unsigned int j=0;pointIdIt != cells.Value()->PointIdsEnd(); pointIdIt++,j++)
			{
				ColorMeshType::PointType pt=0;
				meshes[i]->GetPoint (*pointIdIt, &pt);
				for(unsigned int k=0;k<3;k++)
					pt[k]=pt[k]*(cell_i/meshes[i]->GetNumberOfCells());

				dist +=	avgPoints[j].EuclideanDistanceTo(pt);
				dist_inv +=avgPoints[avgPoints.size()-j-1].EuclideanDistanceTo(pt);
			}
			pointIdIt  =cells.Value()->PointIdsBegin();
			for(unsigned int j=0;pointIdIt != cells.Value()->PointIdsEnd(); pointIdIt++,j++)
			{	
				ColorMeshType::PointType pt=0;
				meshes[i]->GetPoint (*pointIdIt, &pt);
				for(unsigned int k=0;k<3;k++)
				{
					if(dist<dist_inv)
						avgPoints[j][k]+=pt[k]/meshes[i]->GetNumberOfCells();
					else
						avgPoints[avgPoints.size()-j-1][k]+=pt[k]/meshes[i]->GetNumberOfCells();
				}
			}
		}
		//clusterCentroids.push_back(avgPoints);
		cells = meshes[i]->GetCells()->Begin();
		float min_dist = std::numeric_limits<float>::max();
		cell_i=0;
		for(;cells!= meshes[i]->GetCells()->End(); cells++)
		{
			ColorMeshType::CellTraits::PointIdIterator  pointIdIt  =cells.Value()->PointIdsBegin();
			double dist=0.0;
			double dist_inv=0.0;
			for(unsigned int j=0;pointIdIt != cells.Value()->PointIdsEnd(); pointIdIt++,j++)
			{
				ColorMeshType::PointType pt=0;
				meshes[i]->GetPoint (*pointIdIt, &pt);

				dist +=	avgPoints[j].EuclideanDistanceTo(pt);
				dist_inv +=avgPoints[avgPoints.size()-j-1].EuclideanDistanceTo(pt);
			}
			if(dist<min_dist)
			{
				min_dist = dist;
				clusterCentroidsIndex[i]= cell_i;
			}
			if(dist_inv<min_dist)
			{
				min_dist = dist_inv;
				clusterCentroidsIndex[i]= cell_i;
			}
		}


	}
	return clusterCentroidsIndex;
}
*/
template <class TColorMesh, class TImage, class THistogramMesh>
void ClusterTools<TColorMesh, TImage, THistogramMesh>::GetPolyDatas(std::vector<std::string> files, std::vector<vtkSmartPointer<vtkPolyData>> polydatas, typename TImage::Pointer image)
{
	for(int i =0; i< files.size(); i++)   
	{
		typename MeshConverterType::Pointer converter = MeshConverterType::New();
		if(  files[i].find(".trk") !=std::string::npos)
		{
			itk::SmartPointer<TrkVTKPolyDataFilter<ImageType>> trkReader  = TrkVTKPolyDataFilter<ImageType>::New();
			trkReader->SetTrkFileName(files[i]);
			trkReader->SetReferenceImage(image);
			trkReader->TrkToVTK();
			polydatas.push_back(trkReader->GetOutputPolyData());

		}
		else
		{
			vtkPolyDataReader *reader = vtkPolyDataReader::New();
			reader->SetFileName ( files[i].c_str());
#if VTK_MAJOR_VERSION > 5
			reader->Update();
			polydatas.push_back(reader->GetOutputPort());
#else
			reader->GetOutput()->Update();
			polydatas.push_back(reader->GetOutput());
#endif
		}

	}
}
template <class TColorMesh, class TImage, class THistogramMesh>
std::vector<typename TColorMesh::Pointer> ClusterTools<TColorMesh, TImage, THistogramMesh>::FixSampleClusters(std::vector<vtkSmartPointer<vtkPolyData>> polydatas, int numberOfPoints)
{

	std::vector<typename ColorMeshType::Pointer> meshes;
	for (unsigned int i=0;i<polydatas.size(); i++)
	{
		vtkSmartPointer<vtkSplineFilter> spline = vtkSmartPointer<vtkSplineFilter>::New();
#if VTK_MAJOR_VERSION > 5
		spline->SetInputData(polydatas[i]);
#else
		spline->SetInput(polydatas[i]);
#endif
		spline->SetNumberOfSubdivisions(numberOfPoints-1);
		spline->Update();

		typename MeshConverterType::Pointer converter = MeshConverterType::New();
		converter->SetVTKPolyData ( spline->GetOutput() );
		converter->GenerateData2();

		typename ColorMeshType::Pointer mesh =  converter->GetOutput();
		meshes.push_back(mesh);	
	}
	return meshes;
}


template <class TColorMesh, class TImage, class THistogramMesh>
int ClusterTools<TColorMesh, TImage, THistogramMesh>::GetAverageStreamline(typename TColorMesh::Pointer mesh)
{
	//typename ColorMeshType::Pointer mesh =  mesh;

	typename ColorMeshType::CellsContainer::ConstIterator cells = mesh->GetCells()->Begin();
	std::vector<typename ColorMeshType::PointType> avgPoints(cells.Value()->GetNumberOfPoints(),0);
	for(;cells!= mesh->GetCells()->End(); cells++)
	{
		typename ColorMeshType::CellTraits::PointIdIterator  pointIdIt;
		double dist=0.0;
		double dist_inv=0.0;
		int j=0;
		for(pointIdIt  =cells.Value()->PointIdsBegin();pointIdIt != cells.Value()->PointIdsEnd(); pointIdIt++,j++)
		{	
			typename ColorMeshType::PointType pt;
			mesh->GetPoint (*pointIdIt, &pt);
			dist +=	avgPoints[j].EuclideanDistanceTo(pt);
			dist_inv +=avgPoints[avgPoints.size()-j-1].EuclideanDistanceTo(pt);
		}
		j=0;
		for(pointIdIt  =cells.Value()->PointIdsBegin();pointIdIt != cells.Value()->PointIdsEnd(); pointIdIt++,j++)
		{	
			typename ColorMeshType::PointType pt;
			mesh->GetPoint (*pointIdIt, &pt);
			for(int k=0;k<3;k++)
			{
				if (dist <dist_inv)
					avgPoints[j][k]+=pt[k]/mesh->GetNumberOfCells();
				else
					avgPoints[avgPoints.size()-j-1][k]+=pt[k]/mesh->GetNumberOfCells();
			}
		}
	}
	cells = mesh->GetCells()->Begin();

	int cellId;

	float minDistance=std::numeric_limits<int>::max();

	int j=0;
	for(;cells!= mesh->GetCells()->End(); cells++)
	{
		typename ColorMeshType::CellTraits::PointIdIterator  pointIdIt;
		double dist=0.0;
		double dist_inv=0.0;
		for(pointIdIt  =cells.Value()->PointIdsBegin();pointIdIt != cells.Value()->PointIdsEnd(); pointIdIt++,j++)
		{	
			typename ColorMeshType::PointType pt;
			mesh->GetPoint (*pointIdIt, &pt);
			dist +=	avgPoints[j].EuclideanDistanceTo(pt);
			dist_inv +=avgPoints[avgPoints.size()-j-1].EuclideanDistanceTo(pt);
		}
		if(minDistance < std::min(dist, dist_inv))
		{
			cellId=j;	
			minDistance =std::min(dist,dist_inv);
		}
		j+=1;
	}
	return cellId;	
}
template <class TColorMesh, class TImage, class THistogramMesh>
float ClusterTools<TColorMesh, TImage, THistogramMesh>::GetStandardDeviation(typename THistogramMesh::Pointer mesh, int clusterMean)
{	

	typename MembershipFunctionType::Pointer function = MembershipFunctionType::New();	
	function->WithCosine(false);	

	typename HistogramMeshType::CellsContainer::ConstIterator cells =mesh->GetCells()->Begin();
	int numberOfPoints= cells.Value()->GetNumberOfPoints();

	MeasurementVectorType average(numberOfPoints*3);	      
	average.SetCell(mesh, clusterMean) ;
	float distance=0;
	for(int j=0;cells!= mesh->GetCells()->End(); cells++, j++)
	{
		MeasurementVectorType mv(numberOfPoints*3);	      
		mv.SetCell(mesh, j) ;
		distance +=pow( 1.0/function->Evaluate(&mv,&average)-1,2);		
	}
	return 	std::sqrt(distance/(mesh->GetNumberOfCells()-1));

}
template <class TColorMesh, class TImage, class THistogramMesh>
float ClusterTools<TColorMesh, TImage, THistogramMesh>::GetDistance(typename THistogramMesh::Pointer mesh, int clusterMean, int cellId)
{	
	typename MembershipFunctionType::Pointer function = MembershipFunctionType::New();	
	function->WithCosine(false);	

	typename HistogramMeshType::CellsContainer::ConstIterator cells = mesh->GetCells()->Begin();
	int numberOfPoints= cells.Value()->GetNumberOfPoints();

	MeasurementVectorType average(numberOfPoints*3);	      
	average.SetCell(mesh, clusterMean) ;
	
	MeasurementVectorType mv(numberOfPoints*3);	      
	mv.SetCell(mesh, cellId) ;
	float distance =pow( 1.0/function->Evaluate(&mv,&average)-1,2);		
	return distance;
}


#endif
