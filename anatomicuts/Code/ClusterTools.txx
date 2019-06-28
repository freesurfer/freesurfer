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
std::vector<itk::Vector<float>>* ClusterTools<TColorMesh, TImage, THistogramMesh>::GetDirections(DirectionsType dir)
{
	std::vector<itk::Vector<float>>* direcciones = new std::vector<itk::Vector<float>>();
	int possibles[3] = {0,1,-1};
	for(int i=0;i<3;i++)
	{
		for(int k=0;k<3;k++)
		{
			for(int j=0;j<3;j++)
			{
				itk::Vector<float> index;
				index[0] = possibles[i];
				index[1] = possibles[j];
				index[2] = possibles[k];
				int howManyZeros=0;
				if(i==0)
					howManyZeros++;
				if(j==0)
					howManyZeros++;
				if(k==0)
					howManyZeros++;
				if((dir == DirectionsType::ALL && howManyZeros!=3)|| (DirectionsType::STRAIGHT == dir&& howManyZeros== 2) || (DirectionsType::DIAGONAL == dir && howManyZeros>=1 &&howManyZeros!=3)) 
				{
					direcciones->push_back(index);
				}
			}
		}

	}

	return direcciones;
}
template <class TColorMesh, class TImage, class THistogramMesh>
std::vector<typename THistogramMesh::Pointer>* ClusterTools<TColorMesh, TImage, THistogramMesh>::ColorMeshToHistogramMesh(std::vector<typename TColorMesh::Pointer> basicMeshes,typename  TImage::Pointer segmentation, bool removeInterHemispheric)
{
	if( ct == NULL)
	{
		FSENV *fsenv = FSENVgetenv();
		char tmpstr[2000];	
		sprintf(tmpstr, "%s/FreeSurferColorLUT.txt", fsenv->FREESURFER_HOME);
		ct = CTABreadASCII(tmpstr);
	}

	std::vector<typename HistogramMeshType::Pointer>* meshes = new std::vector<typename HistogramMeshType::Pointer>();
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
			meshes->push_back(HistogramMeshType::New());
		}
		else
		{
			meshes->push_back(mesh);
		}
	}
	return meshes;
}

template <class TColorMesh, class TImage, class THistogramMesh>
void ClusterTools<TColorMesh, TImage, THistogramMesh>::SetDirectionalNeighbors(std::vector<typename HistogramMeshType::Pointer>* meshes, typename ImageType::Pointer segmentation, std::vector<itk::Vector<float>> direcciones, bool symmetry)
{
	std::vector<MeasurementVectorType> measurements;
	for(unsigned int i=0;i<meshes->size();i++)
	{	
		int  pointId =0;
		//ColorMeshType::CellTraits::PointIdIterator  pointIdEnd =meshes[i]->GetNumberOfPoints();
		int  pointIdEnd =(*meshes)[i]->GetNumberOfPoints();
		int numPoints = (*meshes)[i]->GetNumberOfPoints();

		for(;pointId != pointIdEnd; pointId++)
		{
			typename HistogramMeshType::PointType pt1=0;
			(*meshes)[i]->GetPoint (pointId, &pt1);
			typename ImageType::IndexType index;	
			if (segmentation->TransformPhysicalPointToIndex(pt1,index))
			{
				HistogramDataType* pointData = new HistogramDataType();
				typename ImageType::PixelType labelOrig = segmentation->GetPixel(index);
				typename ImageType::PixelType label = labelOrig;
				if (symmetry)	
				{
					label= SymmetricLabelId(labelOrig);
				} 

				pointData->push_back(label);


				for(unsigned int k=0;k<direcciones.size();k++)
				{
					typename ImageType::PixelType vecino = labelOrig;
					itk::ContinuousIndex<float,3> continuousIndex = index;
					HistogramPointType point = pt1;
					//std::cout << direcciones[k] << std::endl;
					typename ImageType::IndexType roundedIndex;
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
				(*meshes)[i]->SetPointData(pointId, *pointData);
			}
			else
			{
				std::cout << "ClusterTools::SetDirectionalNeighbors: Point outside image. " << index<<std::endl;
			}
		}
	}

}
template <class TColorMesh, class TImage, class THistogramMesh>
void ClusterTools<TColorMesh, TImage, THistogramMesh>::GetPolyDatas(std::vector<std::string> files, std::vector<vtkSmartPointer<vtkPolyData>>* polydatas, typename TImage::Pointer image)
{
	for(int i =0; i< files.size(); i++)   
	{
		if(  files[i].find(".trk") !=std::string::npos)
		{
			itk::SmartPointer<TrkVTKPolyDataFilter<ImageType>> trkReader  = TrkVTKPolyDataFilter<ImageType>::New();
			trkReader->SetTrkFileName(files[i]);
			trkReader->SetReferenceImage(image);
			trkReader->TrkToVTK();
			polydatas->push_back(trkReader->GetOutputPolyData());

		}
		else
		{
			vtkPolyDataReader *reader = vtkPolyDataReader::New();
			reader->SetFileName ( files[i].c_str());
#if VTK_MAJOR_VERSION > 5
			reader->Update();
#else
			reader->GetOutput()->Update();
#endif
			polydatas->push_back(reader->GetOutput());
		}

	}
}

template <class TColorMesh, class TImage, class THistogramMesh>
void ClusterTools<TColorMesh, TImage, THistogramMesh>::SaveMesh(typename TColorMesh::Pointer mesh ,typename TImage::Pointer image,  std::string outputName, std::string refFiber)
{

	typename VTKConverterType::Pointer vtkConverter =  VTKConverterType::New();
	vtkConverter->SetInput(mesh);
	vtkConverter->Update();
	if(  std::string(outputName).find(".trk") != std::string::npos)
	{
		itk::SmartPointer<TrkVTKPolyDataFilter<ImageType>> trkReader  = TrkVTKPolyDataFilter<ImageType>::New();
		trkReader->SetInput(vtkConverter->GetOutputPolyData());
		trkReader->SetReferenceImage(image);
		trkReader->SetReferenceTrack(	std::string(refFiber));
		trkReader->VTKToTrk(std::string(outputName));

	}
	else
	{
		vtkPolyDataWriter *writerFixed = vtkPolyDataWriter::New();
		writerFixed->SetFileName ( outputName.c_str());
#if VTK_MAJOR_VERSION > 5
		writerFixed->SetInputData(vtkConverter->GetOutputPolyData());
#else
		writerFixed->SetInput(vtkConverter->GetOutputPolyData());
#endif
		writerFixed->SetFileTypeToBinary();
		writerFixed->Update();
	}


}

template <class TColorMesh, class TImage, class THistogramMesh>
std::vector<typename TColorMesh::Pointer>* ClusterTools<TColorMesh, TImage, THistogramMesh>::PolydataToMesh(std::vector<vtkSmartPointer<vtkPolyData>> polydatas)
{
	std::vector<typename ColorMeshType::Pointer>* meshes = new std::vector<typename ColorMeshType::Pointer>();
	for (unsigned int i=0;i<polydatas.size(); i++)
	{
		typename MeshConverterType::Pointer converter = MeshConverterType::New();
		converter->SetVTKPolyData ( polydatas[i] );
		converter->GenerateData2();

		typename ColorMeshType::Pointer mesh =  converter->GetOutput();
		meshes->push_back(mesh);	
	}
	return meshes;
}


template <class TColorMesh, class TImage, class THistogramMesh>
std::vector<typename TColorMesh::Pointer>* ClusterTools<TColorMesh, TImage, THistogramMesh>::FixSampleClusters(std::vector<vtkSmartPointer<vtkPolyData>> polydatas, int numberOfPoints)
{

	std::vector<typename ColorMeshType::Pointer>* meshes = new std::vector<typename ColorMeshType::Pointer>();
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
		meshes->push_back(mesh);	
	}
	return meshes;
}


template <class TColorMesh, class TImage, class THistogramMesh>
int ClusterTools<TColorMesh, TImage, THistogramMesh>::GetAverageStreamline(typename TColorMesh::Pointer mesh)
{
	//typename ColorMeshType::Pointer mesh =  mesh;

	typename ColorMeshType::CellsContainer::ConstIterator cells = mesh->GetCells()->Begin();
	std::vector<typename ColorMeshType::PointType> avgPoints(cells.Value()->GetNumberOfPoints(),0);
	for(int i=0;i<cells.Value()->GetNumberOfPoints();i++)
	{
		typename ColorMeshType::PointType pt;
		pt.Fill(0);
		avgPoints[i]=	pt;
	}
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

	int cellId=-1;

	float minDistance=std::numeric_limits<int>::max();

	for(int k=0;cells!= mesh->GetCells()->End(); cells++, k++)
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
		if(minDistance > std::min(dist, dist_inv))
		{
			cellId=k;	
			minDistance =std::min(dist,dist_inv);
		}
	}
	return cellId;	
}
template <class TColorMesh, class TImage, class THistogramMesh>
float ClusterTools<TColorMesh, TImage, THistogramMesh>::GetStandardDeviation(typename THistogramMesh::Pointer mesh, int clusterMean)
{	

	typename MembershipFunctionType::Pointer function = MembershipFunctionType::New();	
	function->SetMeanEuclidean(true);
	//function->SetLabels(true);
	typename HistogramMeshType::CellsContainer::ConstIterator cells =mesh->GetCells()->Begin();
	int numberOfPoints= cells.Value()->GetNumberOfPoints();

	MeasurementVectorType average(numberOfPoints*3);	      
	average.SetCell(mesh, clusterMean) ;
	float distance=0;
	for(int j=0;cells!= mesh->GetCells()->End(); cells++, j++)
	{
		MeasurementVectorType mv(numberOfPoints*3);	      
		mv.SetCell(mesh, j) ;
		//std::cout << function->Evaluate(&mv, &average) << std::endl;
		distance +=pow( 1.0/function->Evaluate(&mv,&average)-1,2);		
		//distance +=pow( function->Evaluate(&mv,&average),2);		
	}
	return 	std::sqrt(distance/(mesh->GetNumberOfCells()-1));

}
template <class TColorMesh, class TImage, class THistogramMesh>
float ClusterTools<TColorMesh, TImage, THistogramMesh>::GetDistance(typename THistogramMesh::Pointer mesh, int clusterMean, int cellId)
{	
	typename MembershipFunctionType::Pointer function = MembershipFunctionType::New();	
	function->SetMeanEuclidean(true);
	//function->SetLabels(true);
	typename HistogramMeshType::CellsContainer::ConstIterator cells = mesh->GetCells()->Begin();
	int numberOfPoints= cells.Value()->GetNumberOfPoints();

	MeasurementVectorType average(numberOfPoints*3);	      
	average.SetCell(mesh, clusterMean) ;
	
	MeasurementVectorType mv(numberOfPoints*3);	      
	mv.SetCell(mesh, cellId) ;
	float distance =pow( 1.0/function->Evaluate(&mv,&average)-1,2);		
	//float distance =pow( function->Evaluate(&mv,&average),2);		
	return distance;
}


#endif
