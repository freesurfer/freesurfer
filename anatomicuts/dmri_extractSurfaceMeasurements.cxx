/* 
 * Author: Viviana Siless
 * Name: dmri_extractSurfaceMeasurements.cxx
 *
 * Description:
 *
 */

//Libraries
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
// Input Splicing
#include "GetPot.h"

// TRK Reading
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include "itkPolylineCell.h"
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <cmath>
#include "itkArray.h"
#include "itkPolylineCell.h"
#include "TrkVTKPolyDataFilter.txx"
#include "itkImage.h"
#include "PolylineMeshToVTKPolyDataFilter.h"
#include "LabelPerPointVariableLengthVector.h"
#include "EuclideanMembershipFunction.h"
#include "ClusterTools.h"
#include "itkDefaultStaticMeshTraits.h"

// Surface Reading
#include "itkImage.h"
#include <map>
#include "itkDefaultStaticMeshTraits.h"
#include "fsSurface.h"
#include "itkTriangleCell.h"
#include <set>
#include "colortab.h"
#include "fsenv.h"
#include "itkVTKPolyDataWriter.h"
#include "itkSmoothingQuadEdgeMeshFilter.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
	
#include "vtkFillHolesFilter.h" 
#include "vtkPolyDataNormals.h"
#include "vtkCellArray.h"
#include "vtkTriangle.h"
#include "vtkDecimatePro.h"
#include "vtkCleanPolyData.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkTriangleFilter.h"

#include "vtkDelaunay3D.h"
#include "macros.h"
#include "mrisurf.h"
#include "mri.h"
#include "vtkKdTreePointLocator.h"
#include "vtkCurvatures.h"
#include "itkTransformFileReader.h"
#include "itkTransformFactoryBase.h"
#include "itkTransformFactory.h"
#include "itkTransformMeshFilter.h"
#include "itkAffineTransform.h"
#include <itkMatrixOffsetTransformBase.h>
#include "itkTransformFileReader.h"
#include "transform.h"


using namespace std;

// HELPER FUNCTIONS
float calculate_mean(vector<float> n);
float calculate_stde(vector<float> n, float mean);
string makeCSV(string dir, string file, string extension);
vtkIdType which_ID(double n1, double n2, vtkIdType ID1, vtkIdType ID2);
vtkSmartPointer<vtkPolyData> FSToVTK(MRIS* surf);

int main(int narg, char* arg[])
{
	GetPot num1(narg, const_cast<char**>(arg));
	
	// Checking for correct parameters
	if ((num1.size() <= 8) or (num1.search(2, "--help", "-h")))
	{
		cerr << "Usage: " << endl
		     << arg[0] << " -i streamlineFile.trk -sl surfaceFile_lh.orig -tl overlayFile_lh.thickness -cl overlayFile_lh.curv" << endl
		     << "-sr surfaceFile_rh.orig -tr overlayFile_rh.thickness -cr overlayFile_rh.curv -o outputDirectory" << endl 
		     << "-rid reference_image (NOTE: only use reference image when FA is not used" << endl
		     << "-ria reference image for anatomical space (NOTE: when diffusion and anatomical spaces are not registered) " << endl
		     << "-t transformation from diffusion to anatomical space " << endl
		     << "-a annotationFile " << endl
		     << "OPTION: -fa <numFiles> <Filename> FA_file.nii.gz ... <Filename> <fileAddress>" << endl;

		return EXIT_FAILURE;
	}

	// Declaration of Variables for Program to Function
	// TRK file Definitions
	enum {Dimension = 3};
	typedef int PixelType;
	const unsigned int PointDimension = 3;
	typedef vector<int> PointDataType;
	const unsigned int MaxTopologicalDimension = 3;
	typedef double CoordinateType;
	typedef double InterpolationWeightType;
	typedef itk::DefaultStaticMeshTraits<
		PointDataType, PointDimension, MaxTopologicalDimension,
		CoordinateType, InterpolationWeightType, PointDataType > MeshTraits;
	typedef itk::Mesh< PixelType, PointDimension, MeshTraits > HistogramMeshType;
	typedef itk::MatrixOffsetTransformBase<double, 3,3> MatrixOffsetTransformBase_double_3_3;

	typedef itk::Image<float, 3> ImageType;
	
	typedef itk::Mesh< PixelType, PointDimension > ColorMeshType;
	typedef ColorMeshType::PointType PointType;
	typedef ColorMeshType::CellType CellType;
	typedef itk::PolylineCell<CellType> PolylineCellType;
	typedef ColorMeshType::CellAutoPointer CellAutoPointer;
	typedef ClusterTools<ColorMeshType, ImageType, HistogramMeshType> ClusterToolsType;
	
	// Surface Defintions
	typedef float CoordType;
	typedef fs::Surface< CoordType, Dimension> SurfType;
	std::map<long long, int> bundlesIndeces;

	// Input Parsing
	vector<string> TRKFiles;	
	for (string inputName = string(num1.follow("", 2, "-i", "-I")); access(inputName.c_str(), 0) == 0; inputName = string(num1.next("")))
		TRKFiles.push_back(inputName);

	const char *fileCorr =num1.follow("output.csv",2,"-p","-P"); 
	// Left Hemisphere
	const char *surfaceFileL = num1.follow("Left Surface File Not Found", "-sl");
	const char *thickFileL   = num1.follow("Left Thickness File Not Found", "-tl");
	const char *curvFileL    = num1.follow("Left Curvature File Not Found", "-cl");
	
	// Right Hemisphere
	const char *surfaceFileR = num1.follow("Right Surface File Not Found", "-sr");
	const char *thickFileR   = num1.follow("Right Thickness File Not Found", "-tr");
	const char *curvFileR    = num1.follow("Right Curvature File Not Found", "-cr");
	
	const char *outputDir    = num1.follow("Output Directory Not Found", "-o");
	
	const char *refImageDiffusion     = num1.follow("Reference Image Not Found", "-rid");
	const char *refImageSurface = num1.follow("Reference Image Not Found", "-ria");
	const char *annotationFileL= num1.follow("Annotation File Not Found", "-al");
	const char *annotationFileR= num1.follow("Annotation File Not Found", "-ar");
	const char *transformationFile= num1.follow("Transformation matrix not found Not Found", "-t");
	//TRANSFORM* trans  = TransformRead(transformationFile);
	//LTA* lta  = LTAread(transformationFile);
	
	FSENV *fsenv = FSENVgetenv();
	char tmpstr[2000];	
	sprintf(tmpstr, "%s/FreeSurferColorLUT.txt", fsenv->FREESURFER_HOME);
	std::cout << " Color table file " << tmpstr << std::endl; 

	COLOR_TABLE* ct = CTABreadASCII(tmpstr);
	
	// Reading in FA file
	vector<ImageType::Pointer> volumes;
	vector<string> image_fileNames;
	vector<ImageType::Pointer> ref_Image;
	MRI *image;
	
	typedef itk::ImageFileReader<ImageType> ImageReaderType;
	ImageReaderType::Pointer readerS = ImageReaderType::New();
	readerS->SetFileName(refImageDiffusion);
	readerS->Update();
	ref_Image.push_back(readerS->GetOutput());	
	readerS = ImageReaderType::New();
	readerS->SetFileName(refImageSurface);
	readerS->Update();
	ref_Image.push_back(readerS->GetOutput());

	image = MRIread(refImageDiffusion);
	
	int numFiles = num1.follow(0, "-fa");
	bool FA_FOUND = num1.search("-fa");
	num1.next("");
	
	if (FA_FOUND and numFiles > 0) 
	{
		for (int i = 0; i < numFiles; i++)
		{
			image_fileNames.push_back(string(num1.next("")));
			const char *inFile = num1.next("");
			typedef itk::ImageFileReader<ImageType> ImageReaderType;
			ImageReaderType::Pointer readerF = ImageReaderType::New();
			readerF->SetFileName(inFile);
			readerF->Update();
			ImageType::Pointer image  = readerF->GetOutput();
			volumes.push_back(image);	
		}
	} 
	
	//Outputting the Files to Ensure the correct files were input
	cerr << endl;
	for (int i = 0; i < TRKFiles.size(); i++)
	{
		cerr << "TRK File " << i + 1 << ":      " << TRKFiles.at(i) << endl;
		bundlesIndeces[(long long) atoll(TRKFiles.at(i).c_str())]=i;	
		
	}

	cerr << "Left Surface:    " << surfaceFileL << endl << "Left Thickness:  " << thickFileL << endl << "Left Curvature:  " << curvFileL << endl 
	     << "Right Surface:   " << surfaceFileR << endl << "Right Thickness: " << thickFileR << endl << "Right Curvature: " << curvFileR << endl
	     << "Output:          " << outputDir << endl << "Reference Image: " << refImageDiffusion << endl << " Reference Image surface: "<< refImageSurface << endl
		<< " Transformation diffusion to surface: " << transformationFile <<endl  ;
	
	if (FA_FOUND)
	{	
		for (int i = 0; i < image_fileNames.size(); i++)
		{
			cerr << "Image " << i + 1 << ":         " << image_fileNames.at(i) << endl;	
		}
	} 

	// Loading the TRK files into a mesh
	vector<ColorMeshType::Pointer>* meshes;
	vector<vtkSmartPointer<vtkPolyData>> polydatas;
	
	ClusterToolsType::Pointer clusterTools = ClusterToolsType::New();
	clusterTools->GetPolyDatas(TRKFiles, &polydatas, ref_Image.at(0));
	meshes = clusterTools->FixSampleClusters(polydatas,10);
	
	//Loading the surface for each hemisphere and metric
	//Left Curvature
	MRI_SURFACE *surfCL;
        surfCL = MRISread(surfaceFileL);
	
	SurfType::Pointer surfaceCL = SurfType::New();
	surfaceCL->Load(&*surfCL);
	
	surfCL = surfaceCL->GetFSSurface(&*surfCL);
	MRISreadCurvature(surfCL, curvFileL);	
	MRISreadAnnotation(surfCL, annotationFileL);
	surfCL->ct = ct;
	//Left Thickness
	MRI_SURFACE *surfTL;
	surfTL = MRISread(surfaceFileL);

	SurfType::Pointer surfaceTL = SurfType::New();
	surfaceTL->Load(&*surfTL);

	surfTL = surfaceTL->GetFSSurface(&*surfTL);

	MRISreadCurvature(surfTL, thickFileL);

	//Right Curvature
	MRI_SURFACE *surfCR;
	surfCR = MRISread(surfaceFileR);

	SurfType::Pointer surfaceCR = SurfType::New();
	surfaceCR->Load(&*surfCR);

	surfCR = surfaceCR->GetFSSurface(&*surfCR);
	MRISreadCurvature(surfCR, curvFileR);
	MRISreadAnnotation(surfCR, annotationFileR);
	surfCR->ct = ct;
	//Right Thickness
	MRI_SURFACE *surfTR;
	surfTR = MRISread(surfaceFileR);

	SurfType::Pointer surfaceTR = SurfType::New();
	surfaceTR->Load(&*surfTR);

	surfTR = surfaceTR->GetFSSurface(&*surfTR);

	MRISreadCurvature(surfTR, thickFileR);

	// Loading the surface into a KdTree - one for each hemisphere
	// LEFT
	vtkSmartPointer<vtkPolyData> surfVTK_L = FSToVTK(surfCL);
	
	vtkSmartPointer<vtkKdTreePointLocator> surfTreeL = vtkSmartPointer<vtkKdTreePointLocator>::New();
	surfTreeL->SetDataSet(surfVTK_L);
	surfTreeL->BuildLocator();	
	
	// RIGHT
	vtkSmartPointer<vtkPolyData> surfVTK_R = FSToVTK(surfCR);

	vtkSmartPointer<vtkKdTreePointLocator> surfTreeR = vtkSmartPointer<vtkKdTreePointLocator>::New();
	surfTreeR->SetDataSet(surfVTK_R);
	surfTreeR->BuildLocator();	

	// The first and last points in both PointType and an array
	PointType firstPt, lastPt, auxPt;
	firstPt.Fill(0);
	lastPt.Fill(0);
	auxPt.Fill(0);
	double firstPt_array[3];	
	double lastPt_array[3];

	ofstream oFile;
	ofstream averageFile;
	averageFile.open(makeCSV(outputDir, "surfaceMeasures",".csv"));

	// Adds the headers to the files and has option for finding FA values
	averageFile << "streamline,curv.start,curv.end,thickness.start,thickness.end";
	if (FA_FOUND)
	{
		for (int a = 0; a < image_fileNames.size(); a++)
			averageFile << ", mean" << image_fileNames.at(a) << ", stde" << image_fileNames.at(a);
	} 
	averageFile << endl;
	system((std::string("mkdir -p ")+std::string(outputDir)+std::string("/surf")).c_str());
	
	std::ifstream file( fileCorr); 
	std::string value;
	getline ( file, value, ',' ); 
	getline ( file, value, ',' ); 
	std::vector<long long> correspondences ;
	while ( file.good() )
	{
		getline ( file, value, ',' );
		// long long v1 = atoll(value.c_str());
		getline ( file, value, ',' ); 
		long long v2 = atoll(value.c_str());
		correspondences.push_back(v2);		
	} 	
	std::cout << " Files " <<  meshes->size() << std::endl;
	// Cycling through the TRK files
	for(int i = 0; i < meshes->size(); i++)
	{ 
		std::vector<float> values = std::vector<float>(20,0.0);
		// Opening output file with a different name for every TRK File
		//
		oFile.open(makeCSV(std::string(outputDir)+std::string("/surf"), TRKFiles.at(i), ".csv"));

		if (not oFile.is_open()) 
		{
			cerr << "Could not open output file" << endl;
			return -1;
		}

		// Adds the headers to the files and has option for finding FA values
		oFile << "streamline,label.start,label.end,curv.start,curv.end,thickness.start,thickness.end";
		if (FA_FOUND)
		{
			for (int a = 0; a < image_fileNames.size(); a++)
				oFile << ",mean" << image_fileNames.at(a) << ",stde" << image_fileNames.at(a);
			for (int a = 0; a < image_fileNames.size(); a++)
				for (int b =0; b<10; b++)
					oFile << ","<< b << "_"<< image_fileNames.at(a);
		} 
		oFile << ","<<endl;

		// Initialization of a new stream for every TRK files
		
		ColorMeshType::Pointer input = (*meshes)[i] ; //[bundlesIndeces[correspondences[i]]];
		ColorMeshType::CellsContainer::Iterator  inputCellIt = input->GetCells()->Begin();
		
		// Cycling through the streams
		int counter = 1;
		for (; inputCellIt != input->GetCells()->End(); ++inputCellIt, ++counter)
		{
			vector<float> meanFA;
			vector<vector<float>> allFA;
			vector<float> stdeFA;
		
			// If there are image files, then find the mean and stde of FA	
			if (FA_FOUND)
			{	
				
				for (int p = 0; p < volumes.size(); p++)
				{
					allFA.push_back(vector<float>());
					// Creating streamline variable and finding first point
					CellType::PointIdIterator it = inputCellIt.Value()->PointIdsBegin();
					input->GetPoint(*it, &firstPt);
					input->GetPoint(*inputCellIt.Value()->PointIdsEnd(), &lastPt);
					float val = (firstPt[0] -lastPt[0])*(firstPt[1]-lastPt[1])*(firstPt[2]-lastPt[2]);
					
					// Getting the FA value at all points
					vector<float> FA_values;
					ImageType::IndexType index;
				
					// Cycling through the points in the stream
					for (; it != inputCellIt.Value()->PointIdsEnd(); it++)
					{	 
						input->GetPoint(*it, &lastPt);

						// If FA options is used, then add the FA values to the vector	
						if (volumes.at(p)->TransformPhysicalPointToIndex(lastPt, index))
						{
							FA_values.push_back(volumes.at(p)->GetPixel(index));	
							if( val>0)
								allFA[p].push_back(volumes.at(p)->GetPixel(index));	
							else
								allFA[p].insert(allFA[p].begin(),volumes.at(p)->GetPixel(index));	
						}	
					}
			
					// Calculating the MeanFA and stdeFA
					meanFA.push_back(calculate_mean(FA_values));
					stdeFA.push_back(calculate_stde(FA_values, meanFA.at(p)));
				}
			// Otherwise find first and last point of a streamline
			} else {
				CellType::PointIdIterator it = inputCellIt.Value()->PointIdsBegin();
				input->GetPoint(*it, &firstPt);

				for (; it != inputCellIt.Value()->PointIdsEnd(); it++)
					input->GetPoint(*it, &lastPt);
			}

			// Changing the point to an index, then the index to the surface
			ImageType::IndexType first_index, last_index;
       /*  		if( num1.search("-t"))
			{
				//TransformSampleReal2(trans, firstPt[0], firstPt[1], firstPt[2],&auxPt[0],&auxPt[1], &auxPt[2]);
				//LTAworldToWorld(lta, firstPt[0], firstPt[1], firstPt[2],&auxPt[0],&auxPt[1], &auxPt[2]);
		       	auxPt.Fill();
				for(int w=0; w<3;w++)
				{
					auxPt[0]+=lta->xforms[0].m_L(0,w) * firstPt[w] ;  
					auxPt[1]+=lta->xforms[0]->m_L[1][w] * firstPt[w] ;  
					auxPt[2]+=lta->xforms[0]->m_L[2][w] * firstPt[w] ;  
				}
				auxPt[0]+=lta->xforms[0]->m_L[0][3] ;  
				auxPt[1]+=lta->xforms[0]->m_L[1][3] ;  
				auxPt[2]+=lta->xforms[0]->m_L[2][3] ;  */
				/*ref_Image.at(1)->TransformPhysicalPointToIndex(auxPt, first_index);
				LTAworldToWorld(lta, lastPt[0], lastPt[1], lastPt[2],&auxPt[0],&auxPt[1], &auxPt[2]);
				ref_Image.at(1)->TransformPhysicalPointToIndex(auxPt, last_index);
                	}
			else
			*/{
			       	ref_Image.at(0)->TransformPhysicalPointToIndex(firstPt, first_index);
				ref_Image.at(0)->TransformPhysicalPointToIndex(lastPt, last_index);
                	}
			MRIvoxelToSurfaceRAS(image, first_index[0], first_index[1], first_index[2], &firstPt_array[0], &firstPt_array[1], &firstPt_array[2]);
                	MRIvoxelToSurfaceRAS(image, last_index[0], last_index[1], last_index[2], &lastPt_array[0], &lastPt_array[1], &lastPt_array[2]);

			// Finding the vertice number
			double distL, distR;
			vtkIdType Left_ID1  = surfTreeL->FindClosestPointWithinRadius(10, firstPt_array, distL);
			vtkIdType Right_ID1 = surfTreeR->FindClosestPointWithinRadius(10, firstPt_array, distR);		
			vtkIdType ID1 = which_ID(distL, distR, Left_ID1, Right_ID1);

			vtkIdType Left_ID2  = surfTreeL->FindClosestPointWithinRadius(10, lastPt_array, distL);
			vtkIdType Right_ID2 = surfTreeR->FindClosestPointWithinRadius(10, lastPt_array, distR);			
			vtkIdType ID2 = which_ID(distL, distR, Left_ID2, Right_ID2);


			// Outputting values to the file
			std::cout << ID1 << " " <<ID2<< std::endl;
			oFile << "StreamLine" << counter << ", ";
			int structure;
			if (ID1 >=0)
			{
				if (ID1 == Left_ID1)
				{
					CTABfindAnnotation(surfCL->ct , surfCL->vertices[ID1].annotation, &structure);
					std::cout <<"L" << ID1 << " " <<   surfCL->vertices[ID1].annotation << " " <<structure<< std::endl;
				}
				else
				{
					CTABfindAnnotation(surfCR->ct , surfCR->vertices[ID1].annotation, &structure);
					std::cout << "R" << ID1 << " " <<   surfCR->vertices[ID1].annotation << " " <<structure<< std::endl;
				}
				oFile << structure << ",";
			}
			else
			{
				oFile << ",";
			}
			if(ID2>=0)
			{
				if (ID2 == Left_ID2)
				{
					CTABfindAnnotation(surfCL->ct , surfCL->vertices[ID2].annotation, &structure);
					std::cout << "L"<< ID2 << " " <<   surfCL->vertices[ID2].annotation << " " <<structure<< std::endl;
				}
				else
				{
					CTABfindAnnotation(surfCR->ct , surfCR->vertices[ID2].annotation, &structure);
					std::cout <<"R"<< ID2 << " " <<   surfCR->vertices[ID2].annotation << " " <<structure<< std::endl;
				}

				oFile << structure << ",";
			}
			else
			{
				oFile << ",";
			}
			if (ID1 >=0)
			{
				if (ID1 == Left_ID1)
				{
					oFile << surfCL->vertices[ID1].curv << ",";
					values[0]+= surfCL->vertices[ID1].curv ;
				}
				else
				{
					oFile << surfCR->vertices[ID1].curv << ",";
					values[0]+= surfCR->vertices[ID1].curv ;
				}
			}
			else
			{
				oFile << ",";
			}
			if (ID2 >=0)
			{
				if (ID2 == Left_ID2)
				{
					oFile << surfCL->vertices[ID2].curv << ",";
					values[1]+= surfCL->vertices[ID2].curv ;

				}
				else
				{
					oFile << surfCR->vertices[ID2].curv << ",";
					values[1]+= surfCR->vertices[ID2].curv ;
				}
			}
			else
			{
				oFile << ",";
			}	

			if (ID1 >=0)
			{
				if (ID1 == Left_ID1)
				{
					oFile << surfTL->vertices[ID1].curv << ","; 
					values[2]+= surfTL->vertices[ID1].curv ;
				}
				else
				{
					oFile << surfTR->vertices[ID1].curv << ",";
					values[2]+= surfTR->vertices[ID1].curv ;
				}
			}
			else
			{
				oFile << ",";
			}
			if (ID2 >=0)
			{
				if (ID2 == Left_ID2)
				{
					oFile << surfTL->vertices[ID2].curv <<",";
					values[3]+= surfTL->vertices[ID2].curv ;
				}		
				else
				{
					oFile << surfTR->vertices[ID2].curv <<",";
					values[3]+= surfTR->vertices[ID2].curv ;
				}
			}		
			else
			{
				oFile << ",";
			}
			if (FA_FOUND)
			{
				for (int m = 0; m < volumes.size(); m++)
					oFile << meanFA.at(m) << "," << stdeFA.at(m)<<",";
				for (int m = 0; m < volumes.size(); m++)
					for (int b = 0; b < allFA[m].size(); b++)
						oFile <<  allFA[m][b]<<",";
			}

			oFile << endl;
		}

		averageFile <<correspondences[i]<<",";
		for(int m=0; m<4; m++)
			averageFile << values[m]/input->GetNumberOfCells()<< "," ;
		averageFile<< endl;

		oFile.close();
	}
		
	oFile.close();
	averageFile.close();
	return EXIT_SUCCESS;
}



/* Function: makeCSV
 * Input: a directory and a file
 * Returns: a string of a file within the directory
 * Does: Takes in a target directory and returns the directory with the name
 * 	 of the file
 * NOTE: used in conjunction with the creating new CSV files and opening them
 */
string makeCSV(string dir, string file, string extension)
{
	int front = file.find_last_of("/");
	int back  = file.find_last_of(".");	
	if( dir.size() >0)	
		dir.append("/");
	dir.append(file.substr(front + 1, back - front - 1));
	dir.append(extension);

	return dir;
}

/* Function: calculate_mean
 * Input: a vector of all the FA values
 * Return: the mean
 * Does: takes all the values and calculates the mean
 */
float calculate_mean(vector<float> n)
{
	float mean = 0;

	for (int i = 0; i < n.size(); i++)
		mean += n.at(i);

	return mean / n.size();
}

/* Function: calculate_stde
 * Input: a vector of all the FA values and the mean
 * Return: the standard deviation
 * Does: takes all the values and calculates the standard deviation
 */
float calculate_stde(vector<float> n, float mean)
{
	float SD = 0;

	for (int i = 0; i < n.size(); i++)
		SD += pow(n.at(i) - mean, 2);

	return sqrt(SD / n.size());
}

/* Function: which_ID
 * Input: the two distances and the two vertice IDs
 * Return: whichever vertice is closer to the point
 * Does: Compares the two distances and returns the vertice of the shorter distance
 */
vtkIdType which_ID(double n1, double n2, vtkIdType ID1, vtkIdType ID2)
{
	if (n1 < n2)
		return ID1;

	return ID2;	
}

//
// Converts a surface to a VTK
//
vtkSmartPointer<vtkPolyData> FSToVTK(MRIS* surf)
{
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();

	for(int i = 0; i < surf->nvertices; i++)
		points->InsertNextPoint(surf->vertices[i].x, surf->vertices[i].y, surf->vertices[i].z);

	for( int i = 0; i < surf->nfaces; i++ ) {
		vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
		
		for(int j = 0;j < 3; j++)
			triangle->GetPointIds()->SetId(j, surf->faces[i].v[j]);

		triangles->InsertNextCell(triangle);
	}	

	vtkSmartPointer<vtkPolyData> vtkSurface = vtkSmartPointer<vtkPolyData>::New();
	vtkSurface->SetPoints(points);
	vtkSurface->SetPolys(triangles);

	return vtkSurface;
}
