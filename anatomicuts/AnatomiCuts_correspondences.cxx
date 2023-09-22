#include <iostream>
#include "itkImage.h"
#include <map>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_vector_ref.h>
#include "itkDefaultStaticMeshTraits.h"
#include "itkMesh.h"
#include <vnl/vnl_hungarian_algorithm.h>

#if ITK_VERSION_MAJOR >= 5
#include <vcl_legacy_aliases.h>
#else
#include <vcl_limits.h>
#endif
#include <vnl/vnl_matrix.h>
#include "itkImage.h"
#include "itkVector.h"
#include "itkMesh.h"

#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "itkDefaultStaticMeshTraits.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "PolylineMeshToVTKPolyDataFilter.h"
#include "VTKPolyDataToPolylineMeshFilter.h"
//#include "AppendBundleFilter.h"
#include "LabelVariableLengthVector.h"
//#include "FixedVTKSamplingFilter.h"
#include "HierarchicalClusteringPruner.h"
#include "LabelPerPointVariableLengthVector.h"
#include "LabelsEntropyAndIntersectionMembershipFunction.h"
#include <set>
#include "GetPot.h"
#include <string>
#include "vtkSplineFilter.h"
#include "OrientationPlanesFromParcellationFilter.h"


#include "colortab.h"
#include "fsenv.h"


typedef std::vector<int>                  PointDataType;
typedef float PixelType;
const unsigned int PointDimension = 3;
const unsigned int MaxTopologicalDimension = 3;
typedef double CoordinateType;
typedef double InterpolationWeightType;
typedef itk::DefaultStaticMeshTraits<
PointDataType, PointDimension, MaxTopologicalDimension,
	CoordinateType, InterpolationWeightType, PointDataType> MeshTraits;
typedef itk::Mesh< PixelType, PointDimension, MeshTraits > MeshType;
typedef itk::Mesh< PixelType, PointDimension > BasicMeshType;
typedef MeshType::CellType              CellType;
typedef CellType::CellAutoPointer       CellAutoPointer;
typedef PolylineMeshToVTKPolyDataFilter<MeshType> VTKConverterType;
typedef VTKPolyDataToPolylineMeshFilter<BasicMeshType> MeshConverterType;
typedef  itk::Image< double,3> ImageType;
typedef ImageType::IndexType IndexType;
typedef LabelPerPointVariableLengthVector<float, MeshType> MeasurementVectorType;	
typedef LabelsEntropyAndIntersectionMembershipFunction<MeasurementVectorType>  MembershipFunctionType;	


static 	std::string LEFT = "Left";
static	std::string RIGHT = "Right";
static 	std::string LEFT2 = "lh";
static	std::string RIGHT2 = "rh";

COLOR_TABLE *ct=NULL;

int SymmetricLabelId(int id)
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
		//std::cout << sval << std::endl;
	}else if( str.find(LEFT2) != std::string::npos)
	{
		char* hola= (char*)str.replace(str.find(LEFT2),LEFT2.length(), RIGHT2).c_str();
		symId = CTABentryNameToIndex(hola, ct); 

	} else	if( str.find(RIGHT) != std::string::npos)
	{
		char* hola= (char*)str.replace(str.find(RIGHT),RIGHT.length(), LEFT).c_str();
		symId = CTABentryNameToIndex(hola, ct); 
		//std::cout << sval << std::endl;
	}else if( str.find(RIGHT2) != std::string::npos)
	{
		char* hola= (char*)str.replace(str.find(RIGHT2),RIGHT2.length(), LEFT2).c_str();
		symId = CTABentryNameToIndex(hola, ct); 

	}

	/*if( id == symId && id!=0 && id !=5001 && id != 5002)
	  {
	  std::cout << "label " << id << std::endl;
	  }*/
	return  symId;
}
std::vector<MeshType::Pointer> BasicMeshToMesh(std::vector<BasicMeshType::Pointer> basicMeshes, ImageType::Pointer segmentation, float  interHemiThreshold )
{
	if( ct == NULL)
	{
		FSENV *fsenv = FSENVgetenv();
		char tmpstr[2000];	
		sprintf(tmpstr, "%s/FreeSurferColorLUT.txt", fsenv->FREESURFER_HOME);
		ct = CTABreadASCII(tmpstr);
	}

	std::vector<MeshType::Pointer> meshes;
	for(unsigned int i=0;i<basicMeshes.size();i++)
	{
		typedef BasicMeshType::CellsContainer::ConstIterator CellIterator;
		int globalIndex=0;
		int indexCell =0;
		// typedef MeshType::PointIdentifier PointIdentifier;
		typedef MeshType::PointDataContainer PointDataContainerType;
		BasicMeshType::Pointer basicMesh = basicMeshes[i];
		MeshType::Pointer mesh = MeshType::New();
		int out =0;
		for(CellIterator cellIt = basicMesh->GetCells()->Begin(); cellIt!= basicMesh->GetCells()->End(); cellIt++)
		{
			PointDataContainerType::Pointer dataContainer = PointDataContainerType::New();
			int withinIndex=0;
			MeshType::CellAutoPointer line;
			BasicMeshType::CellTraits::PointIdIterator  pointIdIt  = cellIt.Value()->PointIdsBegin();

			BasicMeshType::PointType pt=0; 
			int left=0, right=0;

			if(interHemiThreshold>0)
			{ 
				for(pointIdIt  = cellIt.Value()->PointIdsBegin();pointIdIt != cellIt.Value()->PointIdsEnd();pointIdIt++)
				{
					basicMesh->GetPoint(*pointIdIt, &pt);
					ImageType::IndexType index;
					if (segmentation->TransformPhysicalPointToIndex(pt,index))
					{
						ImageType::PixelType label = segmentation->GetPixel(index);

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
			if (right == 0 || left ==0 || interHemiThreshold<=0 )
			{
				for(pointIdIt  = cellIt.Value()->PointIdsBegin();pointIdIt != cellIt.Value()->PointIdsEnd();pointIdIt++)
				{
					basicMesh->GetPoint(*pointIdIt, &pt);

					line.TakeOwnership ( new itk::PolylineCell<MeshType::CellType> );
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
		if( val> interHemiThreshold)
		{
			meshes.push_back(MeshType::New());
		}
		else
		{
			meshes.push_back(mesh);
		}
	}
	return meshes;
}

std::vector<MeasurementVectorType> SetDirectionalNeighbors(std::vector<MeshType::Pointer> meshes, std::vector<int> clusterCentroidsIndex, ImageType::Pointer segmentation, std::vector<itk::Vector<float>> direcciones, bool symmetry)
{
	//std::cout << " symmetry " << symmetry << std::endl;
	std::vector<MeasurementVectorType> measurements;
	for(unsigned int i=0;i<meshes.size();i++)
	{	
		int  pointId =0;
		//MeshType::CellTraits::PointIdIterator  pointIdEnd =meshes[i]->GetNumberOfPoints();
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
			MeshType::PointType pt1=0;
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
					MeshType::PointType point = pt1;
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
std::vector<int> GetCentroidIndices(std::vector<MeshType::Pointer> meshes)
{
	std::vector<int> clusterCentroidsIndex( meshes.size(),0);
	for(unsigned int i=0; i<meshes.size();i++)	
	{
		std::vector<MeshType::PointType> avgPoints(meshes[i]->GetCells()->Begin().Value()->GetNumberOfPoints(),0);
		MeshType::CellsContainer::ConstIterator cells = meshes[i]->GetCells()->Begin();
		int cell_i=0;
		for(;cells!= meshes[i]->GetCells()->End(); cells++)
		{
			cell_i++;
			MeshType::CellTraits::PointIdIterator  pointIdIt  =cells.Value()->PointIdsBegin();
			double dist=0.0;
			double dist_inv=0.0;
			for(unsigned int j=0;pointIdIt != cells.Value()->PointIdsEnd(); pointIdIt++,j++)
			{
				MeshType::PointType pt=0;
				meshes[i]->GetPoint (*pointIdIt, &pt);
				for(unsigned int k=0;k<3;k++)
					pt[k]=pt[k]*(cell_i/meshes[i]->GetNumberOfCells());

				dist +=	avgPoints[j].EuclideanDistanceTo(pt);
				dist_inv +=avgPoints[avgPoints.size()-j-1].EuclideanDistanceTo(pt);
			}
			pointIdIt  =cells.Value()->PointIdsBegin();
			for(unsigned int j=0;pointIdIt != cells.Value()->PointIdsEnd(); pointIdIt++,j++)
			{	
				MeshType::PointType pt=0;
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
			MeshType::CellTraits::PointIdIterator  pointIdIt  =cells.Value()->PointIdsBegin();
			double dist=0.0;
			double dist_inv=0.0;
			for(unsigned int j=0;pointIdIt != cells.Value()->PointIdsEnd(); pointIdIt++,j++)
			{
				MeshType::PointType pt=0;
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
void GetMeshes(GetPot cl, const char* find1, const char* find2, std::vector<BasicMeshType::Pointer>* meshes, std::vector<std::string>* files)
{
	int maxK =  cl.follow(0,2,find1,find2);

	int k=0;
	while(k<maxK)  
	{

		vtkPolyDataReader *reader = vtkPolyDataReader::New();
		reader->SetFileName ( cl.next("") );
		files->push_back(reader->GetFileName());
#if VTK_MAJOR_VERSION > 5
		reader->Update();
#else
		reader->GetOutput()->Update();
#endif

		vtkSmartPointer<vtkSplineFilter> spline = vtkSmartPointer<vtkSplineFilter>::New();
#if VTK_MAJOR_VERSION > 5
		spline->SetInputConnection(reader->GetOutputPort());
#else
		spline->SetInput(reader->GetOutput());
#endif
		spline->SetNumberOfSubdivisions(9);
		spline->Update();

		MeshConverterType::Pointer converter = MeshConverterType::New();
		converter->SetVTKPolyData ( spline->GetOutput());
		//converter->SetVTKPolyData ( reader->GetOutput() );
		reader->Delete();
		//		converter->Update();
		converter->GenerateData2();



		BasicMeshType::Pointer mesh =  converter->GetOutput();
		meshes->push_back(mesh);	
		k++;

	}
}
std::vector<BasicMeshType::Pointer> FixSampleClusters(std::vector<vtkSmartPointer<vtkPolyData>> polydatas)
{

	std::vector<BasicMeshType::Pointer> meshes;
	for (unsigned int i=0;i<polydatas.size(); i++)
	{
		vtkSmartPointer<vtkSplineFilter> spline = vtkSmartPointer<vtkSplineFilter>::New();
#if VTK_MAJOR_VERSION > 5
		spline->SetInputData(polydatas[i]);
#else
		spline->SetInput(polydatas[i]);
#endif
		spline->SetNumberOfSubdivisions(9);
		spline->Update();

		MeshConverterType::Pointer converter = MeshConverterType::New();
		converter->SetVTKPolyData ( spline->GetOutput() );
		converter->GenerateData2();

		BasicMeshType::Pointer mesh =  converter->GetOutput();
		meshes.push_back(mesh);	
	}
	return meshes;
}
int main(int narg, char*  arg[])
{
	try 
	{
		enum {Dimension =3};
		typedef double                                                        PixelType;
		typedef itk::Image< PixelType,  Dimension >  ImageType;


		GetPot cl(narg, const_cast<char**>(arg));
		GetPot cl2(narg, const_cast<char**>(arg));
		if(cl.size()==1 || cl.search(2,"--help","-h"))
		{
			std::cout<<"Usage: " << std::endl;
			std::cout<< arg[0] << " -s1 parcellation1 -s2 parcellation2 -c numClusters -h1 clusteringPath1  -h2 clusterinPath2 -labels (-euclid for Eucildean) -bb -sym interHemiRatioClusterRemoval -o output"  << std::endl;   
			return -1;
		}
		int numClusters = cl.follow(0,"-c");
		const char *segFile = cl.follow ("", "-s1");
		typedef itk::ImageFileReader<ImageType> ImageReaderType;
		ImageReaderType::Pointer readerS = ImageReaderType::New();
		readerS->SetFileName ( segFile);
		readerS->Update();
		ImageType::Pointer segmentation1  = readerS->GetOutput();



		segFile = cl.follow ("", "-s2");
		float symm =  cl.follow(0.0,"-sym");
		std::cout << "Symmetry " << symm << std::endl;
		bool bb =  cl.search("-bb");
		std::cout << "Baby Mode " << bb << std::endl;
		readerS = ImageReaderType::New();
		readerS->SetFileName ( segFile);
		readerS->Update();
		ImageType::Pointer segmentation2  = readerS->GetOutput();

		std::string hierarchyFilename  = std::string(cl.follow("","-h1"));

		HierarchicalClusteringPruner<BasicMeshType, ImageType>::Pointer hierarchyPruner =  HierarchicalClusteringPruner<BasicMeshType,ImageType>::New();
		hierarchyPruner->SetNumberOfClusters(numClusters);
		hierarchyPruner->SetHierarchyFilename(hierarchyFilename+"/HierarchicalHistory.csv"); 
		hierarchyPruner->SetClustersPath(hierarchyFilename);
		hierarchyPruner->SetExtension(FiberFormat::TRK);
		hierarchyPruner->SetReferenceImage(segmentation1);
		hierarchyPruner->Update();
		std::vector<long long> meshes1Files = hierarchyPruner->GetClustersIds(); 
		std::vector<BasicMeshType::Pointer> basicMeshes1= FixSampleClusters(hierarchyPruner->GetOutputBundles());	
		std::vector<MeshType::Pointer> meshes1= BasicMeshToMesh(basicMeshes1, segmentation1, symm);

		hierarchyFilename  = std::string(cl.follow("","-h2"));

		hierarchyPruner->SetHierarchyFilename(hierarchyFilename+"/HierarchicalHistory.csv"); 
		hierarchyPruner->SetReferenceImage(segmentation2);
		hierarchyPruner->SetClustersPath(hierarchyFilename);
		hierarchyPruner->Update();
		std::vector<BasicMeshType::Pointer> basicMeshes2 = FixSampleClusters(hierarchyPruner->GetOutputBundles());	
		std::vector<MeshType::Pointer> meshes2 = BasicMeshToMesh(basicMeshes2, segmentation2, symm);


		std::vector<long long> meshes2Files = hierarchyPruner->GetClustersIds();
		/*GetMeshes(cl, "-f1" , "-F1", &meshes1,&meshes1Files);
		  std::cout << " meshses 1  size " << meshes1.size() << std::endl;
		  GetMeshes(cl, "-f2" , "-F2", &meshes2, &meshes2Files);
		  std::cout << " meshses 2  size " << meshes2.size() << std::endl;
		  */

		std::vector<int> clusterCentroidsIndex1 = std::vector<int>(meshes1.size(),-1);
		std::vector<int> clusterCentroidsIndex2 =std::vector<int>(meshes2.size(),-1);


		std::map<long long,long long> correspondances;
		vnl_matrix<double> distances(meshes1.size(),meshes2.size());
		if (cl.search(1,"-euclid"))
		{
			std::vector<int> clusterCentroidsIndex1 = GetCentroidIndices(meshes1);
			std::vector<int> clusterCentroidsIndex2 = GetCentroidIndices(meshes2);

			for(unsigned int i=0;i<meshes1.size();i++)
			{	
				float dist_min = std::numeric_limits<float>::max();

				CellAutoPointer cell1;
				meshes1[i]->GetCell(clusterCentroidsIndex1[i], cell1);
				for(unsigned int j=0;j<meshes2.size();j++)
				{
					CellAutoPointer cell2;
					meshes2[j]->GetCell(clusterCentroidsIndex2[j], cell2);

					double dist=0.0;
					double dist_inv=0.0;
					MeshType::CellTraits::PointIdIterator  pointIdIt1  =cell1->PointIdsBegin();
					MeshType::CellTraits::PointIdIterator  pointIdIt2  =cell2->PointIdsBegin();
					MeshType::CellTraits::PointIdIterator  pointIdIt2_inv  =cell2->PointIdsEnd();
					pointIdIt2_inv--;
					for(;pointIdIt1 != cell1->PointIdsEnd(); pointIdIt1++, pointIdIt2++, pointIdIt2_inv--)
					{
						MeshType::PointType pt1=0;
						meshes1[i]->GetPoint (*pointIdIt1, &pt1);
						MeshType::PointType pt2=0;
						meshes2[j]->GetPoint (*pointIdIt2, &pt2);
						MeshType::PointType pt2_inv=0;
						meshes2[j]->GetPoint (*pointIdIt2_inv, &pt2_inv);

						dist +=	pt1.EuclideanDistanceTo(pt2);
						dist_inv += pt1.EuclideanDistanceTo(pt2_inv);
					}
					if(dist<dist_min || dist_inv < dist_min)
					{
						dist_min = std::min(dist, dist_inv);
						correspondances[i]=j;
					}
					distances(i,j)=std::min(dist,dist_inv);

				}
			}	
		}
		else// if (cl.search(1,"-labels"))
		{
			OrientationPlanesFromParcellationFilter<ImageType,ImageType>::Pointer orientationFilter1 = OrientationPlanesFromParcellationFilter<ImageType, ImageType>::New();
			orientationFilter1->SetInput(segmentation1);
			orientationFilter1->SetBabyMode(bb);
			orientationFilter1->Update();
			std::vector<itk::Vector<float>> orientations1;
			orientations1.push_back(orientationFilter1->GetUpDown());
			orientations1.push_back(orientationFilter1->GetFrontBack());
			orientations1.push_back(orientationFilter1->GetLeftRight());

			OrientationPlanesFromParcellationFilter<ImageType,ImageType>::Pointer orientationFilter2 = OrientationPlanesFromParcellationFilter<ImageType, ImageType>::New();
			orientationFilter2->SetInput(segmentation2);
			orientationFilter2->SetBabyMode(bb);
			orientationFilter2->Update();
			std::vector<itk::Vector<float>> orientations2;
			orientations2.push_back(orientationFilter2->GetUpDown());
			orientations2.push_back(orientationFilter2->GetFrontBack());
			orientations2.push_back(orientationFilter2->GetLeftRight());

			//std::cout << orientationFilter1->GetUpDown() << " " << orientationFilter1->GetFrontBack() << " " << orientationFilter1->GetLeftRight() << std::endl;	
			//std::cout << orientationFilter2->GetUpDown() << " " << orientationFilter2->GetFrontBack() << " " << orientationFilter2->GetLeftRight() << std::endl;	


			//typedef ImageType::IndexType IndexType;
			std::vector<itk::Vector<float>> direcciones1;
			std::vector<itk::Vector<float>> direcciones2;
			int possibles[3] = {0,1,-1};
			for(unsigned int i=0;i<3;i++)
			{
				for(unsigned int k=0;k<3;k++)
				{
					for(unsigned int j=0;j<3;j++)
					{
						/*IndexType index1;
						  index1[0] = possibles[i] * orientations1[0];
						  index1[1] = possibles[j] * orientations1[1];
						  index1[2] = possibles[k] * orientations1[2];
						  */
						/*IndexType index2;
						  index2[0] = possibles[i] * orientations2[0];
						  index2[1] = possibles[j] * orientations2[1];
						  index2[2] = possibles[k] * orientations2[2];
						  */	
						itk::Vector<float> dir1;
						itk::Vector<float> dir2;
						for(int w=0;w<3;w++)
						{

							dir1[w] = possibles[i] * orientations1[0][w]+possibles[j] * orientations1[1][w]+possibles[k] * orientations1[2][w];
							if( symm)
							{
								dir2[w] = possibles[i] * orientations2[0][w]+possibles[j] * orientations2[1][w]-possibles[k] * orientations2[2][w];
							}
							else
							{
								dir2[w] = possibles[i] * orientations2[0][w]+possibles[j] * orientations2[1][w]+possibles[k] * orientations2[2][w];
							}
						}
						dir1.Normalize();
						dir2.Normalize();
						int howManyZeros=0;
						if(i==0)
							howManyZeros++;
						if(j==0)
							howManyZeros++;
						if(k==0)
							howManyZeros++;
						if(howManyZeros!=3)
						{
							direcciones1.push_back(dir1);
							direcciones2.push_back(dir2);
							//	std::cout << dir1 << " "<< dir2 << std::endl;
						}
					}
				}
			}
			std::vector<MeasurementVectorType> measurements1 =  SetDirectionalNeighbors(meshes1, clusterCentroidsIndex1,segmentation1, direcciones1, false);
			std::vector<MeasurementVectorType> measurements2 = SetDirectionalNeighbors(meshes2, clusterCentroidsIndex2,segmentation2, direcciones2, symm);
			MembershipFunctionType::Pointer function = MembershipFunctionType::New();		
			function->SetLabels(true);
			//		std::cout << measurements1[100] << std::endl;
			//		std::cout << measurements2[50] << std::endl;

			for(unsigned int i=0;i<measurements1.size();i++)
			{
				float aff_max = std::numeric_limits<float>::min();
				for(unsigned int j=0;j<measurements2.size();j++)
				{
					float aff = function->Evaluate(&measurements1[i], &measurements2[j]);
					//				std::cout << aff << std::endl;
					if(aff > aff_max)
					{
						aff_max = aff;
						correspondances[i]=j;
					}
					//distances(i,j)= (1-aff)*100;
					distances(i,j)= (1/(aff+0.01));

					//distances(i,j)= aff;
				}

			}
		}
		//if(cl.search(1,"-hungarian"))
		{
			//std::cout << " distances " << distances << std::endl;
			vcl_vector<unsigned int> assign = vnl_hungarian_algorithm< double>( distances ); //.GetAssignmentVector();
			//std::cout << " assign size" << assign.size() << std::endl;
			//std::cout << assign[0] << " " << assign[199] << std::endl;
			for( unsigned int i=0;i<assign.size();i++)
				correspondances[i]=assign[i];
		}


		const char *output= cl.follow ("", "-o");

		std::cout << "Correspondances file " << output << std::endl;
		std::ofstream csv_file;
		csv_file.open (output);

		csv_file << "Subject A,Subject B , "<< std::endl;

		//for (std::map<int,std::vector<int>>::iterator it=correspondances.begin(); it!=correspondances.end(); ++it)
		for (std::map<long long,long long>::iterator it=correspondances.begin(); it!=correspondances.end(); ++it)
		{
			csv_file << meshes1Files[it->first] << "," << meshes2Files[it->second] << ","<< distances[it->first][it->second] << std::endl;
		}
		csv_file.close();

	}catch(...)
	{
		std::cout << "Error --> ";
		for(int i=0;i<narg;i++)
		{
			std::cout << arg[i];
		}
		std::cout << std::endl;

	}
	return 0;
}
/*std::cout << "starting hungarian algorithm"<< std::endl;

//subtract smallest entry in row to the row
for(unsigned int i=0;i<distances.rows();i++)
{
vnl_vector<double> row = distances.get_row(i);
row -= row.min_value();
distances.set_row(i, row);
}
//subtract smallest entry in colum to the column
for(unsigned int i=0;i<distances.cols();i++)
{
vnl_vector<double> col = distances.get_column(i);
if( col.min_value() >0 )
{	
col -= col.min_value();
distances.set_column(i, col);
}
}
int covers=0;
do{
std::vector<bool> coveredRows(distances.rows(), false);
std::vector<bool> coveredCols(distances.cols(), false);
covers=0;
//draw minimum amount of lines (in columns and rows) to covered every zero.
for(unsigned int i=0;i<distances.rows();i++)
{
for(unsigned int j=0;j<distances.cols();j++)
{
if(distances(i,j)==0)	
{
int rowZeros=0, colZeros=0;
if(!coveredRows[i]  && !coveredCols[j])
{ 
for(unsigned int k=0;k<distances.cols();k++)
{	
if(distances(i,k)==0 && !coveredCols[k])
rowZeros++;
if(distances(k,j)==0 && !coveredRows[k])
colZeros++;
}
if(rowZeros>colZeros)
coveredRows[i]=true;
else
coveredCols[j]=true;
covers++;
}
}
}
}
//if the covering lines is equal or greater to the matrix size, we finished.
//if not continue
if(covers< distances.rows())
{
//determine  the minimum uncovered entry
double minimum = std::numeric_limits<double>::max();
for(unsigned int i=0;i<distances.rows();i++)
{
for(unsigned int j=0;j<distances.cols();j++)
{	
if(!coveredRows[i] && !coveredCols[j] && (distances(i,j) <minimum ))
minimum= distances(i,j);

}
}

std::cout << " minimum uncovered element "<< minimum << std::endl;
for(unsigned int i=0;i<distances.rows();i++)
{
for(unsigned int j=0;j<distances.cols();j++)
{	
if(coveredCols[j])
	distances(i,j)+=minimum;
if(coveredRows[i])
	distances(i,j)+=minimum;

	}
}

//subtract minimum to every uncoivered row
distances-=minimum;
}
std::cout << "covered lines " << covers << " " << distances.rows() << std::endl;
}while(covers<distances.rows());

std::cout << "getting possible solutions "<< std::endl;
std::map<int, std::set<int>> solution;
for(unsigned int i=0;i<distances.rows();i++)
{
	for(unsigned int j=0;j<distances.cols();j++)
	{
		if(distances(i,j)<=0)
		{
			if(solution.count(i)==0)
				solution[i]= std::set<int>();
			solution[i].insert(j);			
		}
	}
}
std::cout << "removing restrictions "<< std::endl;	
bool foundSolution = false;
bool changed = true;
//remove options when other row needs it
while(changed)
{
	changed=false;
	for(unsigned int i=0;i<distances.rows();i++)
	{
		if(solution[i].size()==1)
		{
			correspondances[i]=*solution[i].begin();
			for(int j=0;j<distances.cols();j++)
			{	if(i!=j && solution[j].count(*solution[i].begin())>0)
				{	
					solution[j].erase(*solution[i].begin());
					changed=true;
				}
			}
		}	
	}
}
//delete options
std::cout << "getting final solution" <<std::endl;
for(unsigned int i=0;i<distances.rows();i++)
{
	auto it = solution[i].begin();
	correspondances[i]=*it;
	if(solution[i].size()>1)
	{
		for(unsigned int j=0;j<distances.cols();j++)
		{	if(i!=j && solution[j].count(*it)>0)
			{	
				solution[j].erase(*it);
			}
		}
	}
}

std::cout << "finish hungarian algorithm "<< std::endl;
*/
