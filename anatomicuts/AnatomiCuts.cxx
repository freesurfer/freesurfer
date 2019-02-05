#include <iostream>
#include <string>         
#include "itkImageFileReader.h"
#include "GetPot.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImage.h"
#include <map>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_vector_ref.h>
#include <queue>
#include <sparse/spMatrix.h>
#include "itkMinimumMaximumImageCalculator.h"
#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "itkDefaultStaticMeshTraits.h"
#include "itkMesh.h"
#include "itkPolylineCell.h"
#include <type_traits>
#include <typeinfo>
#define PI 3.14159265
#include "itkRigid3DTransform.h"
#include "NormalizedCutsFilter.h"
#include "TrkVTKPolyDataFilter.txx"
#include "PolylineMeshToVTKPolyDataFilter.h"
#include "itkImageDuplicator.h"
#include "itkNeighborhoodIterator.h"
#include <time.h>
#include "vtkSplineFilter.h"
#include "VTKPolyDataToPolylineMeshFilter.h"
#include "EuclideanMembershipFunction.h"
#include "HausdorffMembershipFunction.h"
#include "LabelPerPointMembershipFunction.h"
#include "LabelsHistogramMembershipFunction.h"
#include "LabelsEntropyAndIntersectionMembershipFunction.h"
#include "LabelsPointToPointMembershipFunction.h"
#include "LabelPerPointVariableLengthVector.h"
#include "vtkDirectory.h" 

typedef struct {
	unsigned char r;
	unsigned char g;
	unsigned char b;
} color_triplet2;

typedef std::pair<int,int> Node;
typedef std::vector<Node> Adjacencies;
typedef std::vector<Node> Graph;

int main(int narg, char* arg[])
{
	enum {Dimension =3};
	typedef int                                                        PixelType;
	typedef itk::Image< PixelType,Dimension> ImageType;
	typedef ImageType::IndexType 			IndexType;
	// typedef itk::MinimumMaximumImageCalculator<ImageType>   MinMaxCalculatorType;

	GetPot cl(narg, const_cast<char**>(arg));
	if(cl.size()==1 || cl.search(2,"--help","-h"))
	{
		std::cout<<"Usage: " << std::endl;
		std::cout<< arg[0] << " -s segmentationFile -f fiber.vtk -c #clusters -n #points  -e #fibers for eigen  -o outputFolder -d [s:straight d:diagonal a:all o:none] "  << std::endl;
		return -1;
	}
	
	const char *segFile = cl.follow ("", "-s");
	const char *fiberFile = cl.follow ("", "-f");
	const char *outputFolder = cl.follow ("", "-o");
	const char *neighbors = cl.follow ("a", "-d");
	int numberOfClusters = cl.follow(200,"-c");
	int numberOfPoints = cl.follow(10, "-n");
	int numberOfFibers = cl.follow(500, "-e");
	vtkDirectory::MakeDirectory(outputFolder);
	std::vector<std::string> labels;
	std::vector<std::pair<std::string,std::string>> clusterIdHierarchy;
	ImageType::Pointer segmentation; 
	std::cout << "number of points  "<< numberOfPoints << std::endl;
	const unsigned int PointDimension = 3;
	{


		typedef std::vector<int>                  PointDataType;
		const unsigned int MaxTopologicalDimension = 3;
		typedef double CoordinateType;
		typedef double InterpolationWeightType;
		typedef itk::DefaultStaticMeshTraits<
			PointDataType, PointDimension, MaxTopologicalDimension,
			CoordinateType, InterpolationWeightType, PointDataType > MeshTraits;
		typedef itk::Mesh< PixelType, PointDimension, MeshTraits > MeshType;

		std::vector<IndexType> direcciones;
		typedef itk::Rigid3DTransform<double>  TransformType;
		std::vector<TransformType::Pointer> rotaciones;
		//bool normal=false;
		std::cout << "neighbors" << neighbors << std::endl;
		int possibles[3] = {0,1,-1};
		for(int i=0;i<3;i++)
		{
			for(int k=0;k<3;k++)
			{
				for(int j=0;j<3;j++)
				{
					IndexType index;
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
					if((neighbors[0] == 'a' && howManyZeros!=3)|| (neighbors[0] == 's' && howManyZeros== 2) || (neighbors[0] == 'd' && howManyZeros>=1 &&howManyZeros!=3)) // ||neighbors[0]== '1' && direcciones.size()==0  )
					{
						direcciones.push_back(index);
					//	std::cout << index << std::endl;
					}
				}
			}

		}
		std::cout << "#directions"<< direcciones.size() << std::endl;

		MeshType::Pointer mesh = MeshType::New();
		mesh->SetCellsAllocationMethod( MeshType::CellsAllocatedDynamicallyCellByCell );			
		//Copy meshes 
		{
			typedef itk::ImageFileReader<ImageType> ImageReaderType;
			ImageReaderType::Pointer reader = ImageReaderType::New();
			reader->SetFileName ( segFile);
			reader->Update();
			segmentation  = reader->GetOutput();
	
			//write a filter
			typedef itk::ImageDuplicator< ImageType > DuplicatorType;
			DuplicatorType::Pointer duplicator = DuplicatorType::New();
			duplicator->SetInputImage(segmentation);
			duplicator->Update();
			ImageType::Pointer clonedImage = duplicator->GetOutput();


			ImageType::SizeType radius;
			radius.Fill(1); 
			int leftWrongPixels=1;
			int iteration=0;
			while(leftWrongPixels>0 && iteration<10)
			{ 
				//std::cout << " iteration " << iteration <<" leftWrongPixels "<<  leftWrongPixels << std::endl;
				leftWrongPixels=0;
				itk::NeighborhoodIterator<ImageType> iter(radius, segmentation,segmentation->GetLargestPossibleRegion());

				for(;!iter.IsAtEnd();++iter)
				{
					if(iter.GetCenterPixel()==0 || iter.GetCenterPixel() ==5001 || iter.GetCenterPixel() ==5002 ||iter.GetCenterPixel() ==41 || iter.GetCenterPixel() ==2  )
					{
						int newPixel = -1;
						for(unsigned int i = 0; i < iter.Size(); i++)
						{
							bool IsInBounds;
							int pixel = iter.GetPixel(i, IsInBounds);
							if(pixel != 0 && pixel!= 5001 &&pixel!=5002)
							{
								newPixel = pixel;
							}
						}
						if(newPixel != -1)
						{
							clonedImage->SetPixel(iter.GetIndex(), newPixel);
						}
						else
						{
							leftWrongPixels++; 
						}
					}
				}

				DuplicatorType::Pointer duplicator = DuplicatorType::New();
				duplicator->SetInputImage(clonedImage);
				duplicator->Update();
				segmentation= duplicator->GetOutput();

				iteration++;
			}
			//end write a filter

			typedef itk::Mesh< PixelType, PointDimension > MeshBasicType;

			typedef VTKPolyDataToPolylineMeshFilter<MeshBasicType> MeshConverterType;
			MeshConverterType::Pointer converter = MeshConverterType::New();

			if(  std::string(fiberFile).find(std::string(".trk")) !=std::string::npos)
			{
				std::cout << "reading trk : " << fiberFile <<" " << std::string(fiberFile).find(".trk")<< std::endl;
				itk::SmartPointer<TrkVTKPolyDataFilter<ImageType>> trkReader  = TrkVTKPolyDataFilter<ImageType>::New();
				trkReader->SetTrkFileName(fiberFile);
				//trkReader->SetReferenceImage(segmentation);
				trkReader->TrkToVTK();
			
				vtkSmartPointer<vtkSplineFilter> spline = vtkSmartPointer<vtkSplineFilter>::New();
				spline->SetInput(trkReader->GetOutputPolyData());
				spline->SetNumberOfSubdivisions(numberOfPoints);
				spline->Update();
				converter->SetVTKPolyData ( spline->GetOutput());

				//converter->SetVTKPolyData ( trkReader->GetOutputPolyData() );
			}
			else
			{
				vtkSmartPointer<vtkPolyDataReader> vtkReader = vtkPolyDataReader::New();
				std::cout <<" reading file: "<<  fiberFile << std::endl;
				vtkReader->SetFileName ( fiberFile);
				vtkReader->Update();
				
				vtkSmartPointer<vtkSplineFilter> spline = vtkSmartPointer<vtkSplineFilter>::New();
				spline->SetInput(vtkReader->GetOutput());
				spline->SetNumberOfSubdivisions(numberOfPoints);
				spline->Update();
				converter->SetVTKPolyData ( spline->GetOutput());
				//converter->SetVTKPolyData ( vtkReader->GetOutput() );
			}	

			//const unsigned int PointDimension = 3;

			converter->Update();

			MeshBasicType::Pointer basicMesh = converter->GetOutput();
			//std::cout << "number of fibers" <<  basicMesh->GetNumberOfCells()<< std::endl;
			/*typedef itk::FixedVTKSamplingFilter<MeshBasicType,MeshBasicType> SamplingFilterType;

			SamplingFilterType::Pointer samplingFilter =  SamplingFilterType::New();
			samplingFilter->SetInput(basicMesh);
			//std::cout <<  "sampling " << std::endl;
			samplingFilter->SetSampling(numberOfPoints);
			basicMesh->GetCells()->Begin();	
			samplingFilter->Update();
			basicMesh = samplingFilter->GetOutput();
			*/
			//std::cout << " finish sampling "<< std::endl;

			typedef MeshBasicType::CellsContainer::ConstIterator CellIterator;
			int globalIndex=0;
			int indexCell =0;
			// typedef MeshType::PointIdentifier PointIdentifier;
			typedef MeshType::PointDataContainer PointDataContainerType;
			//int outsidePoints = 0;
			int numCellsPase =0;

			for(CellIterator cellIt = basicMesh->GetCells()->Begin(); cellIt!= basicMesh->GetCells()->End(); cellIt++)
			{
				PointDataContainerType::Pointer dataContainer = PointDataContainerType::New();
				int withinIndex=0;
				MeshType::CellAutoPointer line;
				MeshBasicType::CellTraits::PointIdIterator  pointIdIt  = cellIt.Value()->PointIdsBegin();

				MeshBasicType::PointType ptMean;
				ptMean.Fill(0.0);
  
				MeshBasicType::PointType pt;  
				int jota=0;
				numCellsPase++;
				for(pointIdIt  = cellIt.Value()->PointIdsBegin();pointIdIt != cellIt.Value()->PointIdsEnd()&&jota<numberOfPoints;pointIdIt++)
				{
					jota++;
					basicMesh->GetPoint(*pointIdIt, &pt);
					
					if(cl.search(1,"-meanAndCov"))
					{
						for(int yy=0;yy<3;yy++)
							ptMean[yy]+=pt[yy];	
					}	
					else
					{	
						line.TakeOwnership ( new itk::PolylineCell<MeshType::CellType> );
						mesh->SetPoint (globalIndex, pt);
						line->SetPointId (withinIndex, globalIndex);

						IndexType index;	
						segmentation->TransformPhysicalPointToIndex(pt,index);				

						PointDataType* pointData = new PointDataType();
						PixelType label = segmentation->GetPixel(index);
						pointData->push_back(label);

						for(unsigned int i=1;i<direcciones.size();i++)
						{
							PixelType vecino = label;
							IndexType ind = index;
							MeshBasicType::PointType point = pt;
							while(vecino == label)
							{
								for(int j=0; j<3;j++)
									ind[j] += direcciones[i][j];
								if(!segmentation->GetLargestPossibleRegion().IsInside(ind))
									break;
								vecino = segmentation->GetPixel(ind);
							}
							if(vecino!=0)	
								pointData->push_back(vecino);
							else
								pointData->push_back(label);

						}
						pointData->push_back(indexCell+9000);	
						mesh->SetPointData(globalIndex, *pointData);
						withinIndex++;
						globalIndex++;
					}
				}
				
				if(cl.search(1,"-meanAndCovInvert") ||cl.search(1,"-meanAndCovGaussian")  )
				{
					for(int yy=0;yy<3;yy++)
						ptMean[yy]/=numberOfPoints;
					line.TakeOwnership ( new itk::PolylineCell<MeshType::CellType> );
					mesh->SetPoint (globalIndex, ptMean);
					line->SetPointId (withinIndex, globalIndex);
					withinIndex++;
					globalIndex++;

					MeshBasicType::PointType ptx ;	
					MeshBasicType::PointType ptyz ;	
					ptx.Fill(0.0);
					ptyz.Fill(0.0);
					
					for(pointIdIt  = cellIt.Value()->PointIdsBegin();pointIdIt != cellIt.Value()->PointIdsEnd();pointIdIt++)
					{
						basicMesh->GetPoint(*pointIdIt, &pt);
					
						MeshBasicType::PointType pt2  = pt -ptMean;	
						ptx[0] += pt2[0]*pt2[0];
						ptx[1] += pt2[0]*pt2[1];
						ptx[2] += pt2[0]*pt2[2];	
						
						ptyz[0] += pt2[1]*pt2[1];
						ptyz[1] += pt2[1]*pt2[2];
						ptyz[2] += pt2[2]*pt2[2];	
						
					}
					
					for(int yy=0;yy<3;yy++)
					{
						ptx[yy]/=(numberOfPoints-1);
						ptyz[yy]/=(numberOfPoints-1);
					}	
					mesh->SetPoint (globalIndex, ptx);
					line->SetPointId (withinIndex, globalIndex);
					withinIndex++;
					globalIndex++;
					mesh->SetPoint (globalIndex, ptyz);
					line->SetPointId (withinIndex, globalIndex);
				}
				else
				{
					//to force dist(x,y)=0 <-> x=y
					mesh->SetPoint (globalIndex, pt);
					line->SetPointId (withinIndex, globalIndex);
					PointDataType* pointData = new PointDataType();
					//for(int i=0;i<direcciones.size();i++)
					//	pointData->push_back(-indexCell-1);	
					mesh->SetPointData(globalIndex, *pointData);
				}
				globalIndex++;
				//end to force dist(x,y)=0 <-> x=y
				
				mesh->SetCell (indexCell, line);
				indexCell++;
			}

			//std::cout << " points mesh " << mesh->GetNumberOfPoints()  << "  basic mesh " << basicMesh->GetNumberOfPoints() <<std::endl;
			//std::cout << " number of fibers " << mesh->GetNumberOfCells()  << std::endl;
		}
		typedef LabelPerPointVariableLengthVector<float, MeshType> MeasurementVectorType;	
		int numMem = numberOfFibers;

		typedef LabelPerPointMembershipFunction<MeasurementVectorType> MembershipFunctionBaseType;	
		std::vector<MembershipFunctionBaseType::Pointer> functionList;

		for(int i =0;i<numMem;i++)
		{
			if(cl.search(1,"-labelsPTPNN") || cl.search(1,"-labelsPTP"))
			{	typedef LabelsPointToPointMembershipFunction<MeasurementVectorType> MembershipFunctionType;
				MembershipFunctionType::Pointer function2 = MembershipFunctionType::New();	
				if( cl.search(1,"-labelsPTP"))
					function2->SetLabelsCount(direcciones.size()+1); // all directions plus the label itself
				else
					function2->SetLabelsCount(1);
				MembershipFunctionBaseType::Pointer function= static_cast<MembershipFunctionBaseType::Pointer>(function2.GetPointer());
				functionList.push_back(function);
			}else if(cl.search(1,"-euclid"))
			{
				//std::cout << "-euclid" << std::endl;
				typedef EuclideanMembershipFunction<MeasurementVectorType> MembershipFunctionType;
				MembershipFunctionType::Pointer function2 = MembershipFunctionType::New();	
				function2->WithCosine(false);	
				MembershipFunctionBaseType::Pointer function= static_cast<MembershipFunctionBaseType::Pointer>(function2.GetPointer());
				functionList.push_back(function);
			}else if(cl.search(1,"-hausdorff"))
			{
				typedef HausdorffMembershipFunction<MeasurementVectorType> MembershipFunctionType;
				MembershipFunctionType::Pointer function2 = MembershipFunctionType::New();	
				function2->WithCosine(false);	
				MembershipFunctionBaseType::Pointer function= static_cast<MembershipFunctionBaseType::Pointer>(function2.GetPointer());
				functionList.push_back(function);
			/*}else if(cl.search(1,"-gaussian"))
			{
				typedef GaussianKernelMembershipFunction<MeasurementVectorType> MembershipFunctionType;
				MembershipFunctionType::Pointer function2 = MembershipFunctionType::New();	
				function2->WithCosine(false);	
				MembershipFunctionBaseType::Pointer function= static_cast<MembershipFunctionBaseType::Pointer>(function2.GetPointer());
				functionList.push_back(function);
			*/
			}else
			{
				typedef LabelsEntropyAndIntersectionMembershipFunction<MeasurementVectorType>  MembershipFunctionType;	
				MembershipFunctionType::Pointer function2 = MembershipFunctionType::New();		
				function2->SetIntersection(cl.search(1,"-intersection"));
				function2->SetEntropy(cl.search(1,"-entropy"));
				function2->SetLabels(cl.search(1,"-labels2") || cl.search(1,"-labels"));
				function2->SetEuclidean(cl.search(1,"-leuclid"));
				function2->SetLabelsAndEuclid(cl.search(1,"-labelsAndEuclid"));
				function2->SetDice(cl.search(1,"-dice"));
				function2->SetKulczynskis(cl.search(1,"-kulczynskis"));

				function2->SetJensenShannon(cl.search(1,"-jensenShannon"));
				function2->SetRuzicka(cl.search(1,"-ruzicka"));
				function2->SetMeanEuclidean(cl.search(1,"-meanEuclidean"));
				function2->SetGaussian(cl.search(1,"-gaussian "));
				function2->SetMeanClosestPointInvert(cl.search(1,"-meanClosestPointInvert"));
				function2->SetMeanClosestPointGaussian(cl.search(1,"-meanClosestPointGaussian"));
				function2->SetMeanAndCovGaussian(cl.search(1,"-meanAndCovGaussian"));
				function2->SetMeanAndCovInvert(cl.search(1,"-meanAndCovInvert"));

				MembershipFunctionBaseType::Pointer function= static_cast<MembershipFunctionBaseType::Pointer>(function2.GetPointer());
				functionList.push_back(function);
			}

		}
		//std::cout << " NormalizedCuts " << std::endl;
		//clock_t t = clock();
		time_t timer1, timer2;
  		double seconds;

  		time(&timer1);  /* get current time; same as: timer = time(NULL)  */


		typedef NormalizedCutsFilter<MeshType, MembershipFunctionBaseType> NormalizeCutsType;
		NormalizeCutsType::Pointer normalizeCuts = NormalizeCutsType::New();
		normalizeCuts->SetNumberOfClusters(numberOfClusters);
		normalizeCuts->SetMembershipFunctionVector(&functionList);
		normalizeCuts->SetNumberOfFibersForEigenDecomposition(numberOfFibers);
		normalizeCuts->SetInput(mesh);
		normalizeCuts->Update();

		labels  = normalizeCuts->GetLabels();
		clusterIdHierarchy = normalizeCuts->GetClusterIdHierarchy();
		time(&timer2);  /* get current time; same as: timer = time(NULL)  */
		seconds = difftime(timer1,timer2);
		//t = clock() - t;
		std::cout << "Execution time: " << seconds/60.0<< " mins"<< std::endl;

	}
	typedef itk::Mesh< PixelType, PointDimension > MeshBasicType;

	typedef VTKPolyDataToPolylineMeshFilter<MeshBasicType> MeshConverterType;
	MeshConverterType::Pointer converter = MeshConverterType::New();
	if(  std::string(fiberFile).find(std::string(".trk")) !=std::string::npos)
	{
		itk::SmartPointer<TrkVTKPolyDataFilter<ImageType>> trkReader  = TrkVTKPolyDataFilter<ImageType>::New();
		trkReader->SetTrkFileName(fiberFile);
		//trkReader->SetReferenceImage(segmentation);
		trkReader->TrkToVTK();
		converter->SetVTKPolyData ( trkReader->GetOutputPolyData() );
	}
	else
	{
		vtkSmartPointer<vtkPolyDataReader> vtkReader = vtkPolyDataReader::New();
		//std::cout <<" reading file: "<<  fiberFile << std::endl;
		vtkReader->SetFileName ( fiberFile);
		vtkReader->Update();
		converter->SetVTKPolyData ( vtkReader->GetOutput() );
	}	
	converter->Update();

	MeshBasicType::Pointer originalMesh= converter->GetOutput();

	std::map<std::string,MeshBasicType::Pointer> newMeshes;
	std::map<std::string,MeshBasicType::PointIdentifier> pointIndices;
	std::map<std::string,MeshBasicType::CellIdentifier> cellIndices;
	//int noClusterFibers=0;

	// typedef itk::PolylineCell<MeshBasicType::CellType>                      PolylineCellType;
	for (unsigned int i=0;i<labels.size();i++)
	{
		MeshBasicType::CellAutoPointer line;
		line.TakeOwnership (new itk::PolylineCell<MeshBasicType::CellType> );
		int k=0;
		MeshBasicType::CellAutoPointer cell;
		originalMesh->GetCell(i, cell);
		MeshBasicType::CellPixelType pixelType=0;
		originalMesh->GetCellData(i, &pixelType);

		if(newMeshes.count(labels[i])==0)
		{
			MeshBasicType::Pointer om = MeshBasicType::New();
			om->SetCellsAllocationMethod( MeshBasicType::CellsAllocatedDynamicallyCellByCell );

			newMeshes[labels[i]]=om;
			pointIndices[labels[i]]=0;
			cellIndices[labels[i]]=0;
		}				


		for (MeshBasicType::CellType::PointIdIterator it = cell->PointIdsBegin(); it!=cell->PointIdsEnd(); it++)
		{
			MeshBasicType::PointType pt;
			pt.Fill (0.0);
			originalMesh->GetPoint (*it, &pt);

			newMeshes[labels[i]]->SetPoint (pointIndices[labels[i]], pt);
			line->SetPointId (k,pointIndices[labels[i]]);

			k++;
			pointIndices[labels[i]]++;

		}
		newMeshes[labels[i]]->SetCell (cellIndices[labels[i]], line);
		newMeshes[labels[i]]->SetCellData (cellIndices[labels[i]], pixelType);


		cellIndices[labels[i]]++;
		//if( labels[i] >=  numberOfClusters)
		//{
		//	noClusterFibers++;
		//}
	}
	//std::cout << " Fibers without clusters  " << noClusterFibers <<std::endl;

	int i;
	int red, green, blue;
	color_triplet2 table[256] = {{0,0,0}}; /* Initialize to all black */

	const color_triplet2 black = {0,0,0};
	const color_triplet2 white = {255,255,255};

	table[0]   = white;
	table[255] = black;


	i = 20; /* first 20 and last 20 are reserved */
	for (red = 0; red <= 255; red+= 51) {/* the six values of red */
		for (green = 0; green <= 255; green += 51) {
			for (blue = 0; blue <= 255; blue+= 51) {
				table[i].r = red;
				table[i].g = green;
				table[i].b = blue;
				++i;
			}
		}
	}
	table[0]   = white; 
	table[255] = black;

	i=0;
	for (std::map<std::string,MeshBasicType::Pointer>::iterator it=newMeshes.begin(); it!=newMeshes.end(); ++it)
	{ 
		int index =((int)47.*((i%newMeshes.size())%(150)))%197+5;
		index = (int)(13*(i%newMeshes.size()))%150+65;
		//std::cout << "index " << index << std::endl;
		unsigned char color[3] = { table[index].r, table[index].g,table[index].b};
		//std::cout << color[0] << " " <<color[1] << " "<<color[2] <<std::endl;
		MeshBasicType::Pointer mesh=it->second; 
		if( mesh->GetNumberOfPoints() >0 )
		{
			typedef PolylineMeshToVTKPolyDataFilter<MeshBasicType> VTKConverterType;

			VTKConverterType::Pointer vtkConverter =  VTKConverterType::New();
			vtkConverter->SetColor(color);
			vtkConverter->SetInput(mesh);
			vtkConverter->Update();
			char meshName2[100];
			sprintf(meshName2, "%s/%s.vtk",outputFolder, it->first.c_str());

			vtkSmartPointer<vtkPolyDataWriter> writerFixed = vtkPolyDataWriter::New();
			writerFixed->SetFileName ( meshName2);
			writerFixed->SetInput(vtkConverter->GetOutputPolyData());
			writerFixed->SetFileTypeToBinary();
			writerFixed->Update();
			//std::cout << " saving file " << meshName2 << std::endl;
	
			itk::SmartPointer<TrkVTKPolyDataFilter<ImageType>> trkReader  = TrkVTKPolyDataFilter<ImageType>::New();
			trkReader->SetInput(vtkConverter->GetOutputPolyData());
			//trkReader->SetReferenceImage(segmentation);
			trkReader->SetReferenceTrack(fiberFile);
			trkReader->SetColor(color);
			sprintf(meshName2, "%s/%s.trk",outputFolder, it->first.c_str());
			trkReader->VTKToTrk(meshName2);

		}
		i++;
	}


	char csv_filename[100];
	sprintf(csv_filename, "%s/HierarchicalHistory.csv",outputFolder);

	std::ofstream csv_file(csv_filename);

	csv_file << "Parent, Child," <<std::endl;
	for (std::vector<std::pair<std::string,std::string>>::iterator it=clusterIdHierarchy.begin(); it!=clusterIdHierarchy.end(); ++it)
		csv_file << (*it).first << "," << (*it).second << ","<<std::endl;

	csv_file.close();
}
