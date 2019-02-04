#include "itkDecisionRule.h"
#include "itkVector.h"
#include "itkListSample.h"
#include "itkKdTree.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkKdTreeBasedKmeansEstimator.h"
#include "itkMinimumDecisionRule.h"
#include "itkEuclideanDistanceMetric.h"
#include "itkDistanceToCentroidMembershipFunction.h"
#include "itkSampleClassifierFilter.h"
#include "itkNormalVariateGenerator.h"

#include "vtkVersion.h"
#include "vtkActor.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkSmartPointer.h"
#include "vtkVertexGlyphFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVariableLengthVector.h"

#include <iostream>
#include "itkImage.h"
#include <map>
#include "itkDefaultStaticMeshTraits.h"
#include "fsSurface.h"
#include "itkTriangleCell.h"
#include <set>
#include "GetPot.h"
#include <string>
#include "colortab.h"
#include "fsenv.h"
#include "fsSurfaceOptimizationFilter.h"
#include "itkVTKPolyDataWriter.h"
#include <vnl/vnl_cross.h>
#include <cmath>

extern "C" 
{
	#include "mrisurf.h"
}

int main(int narg, char * arg[])
{
	constexpr unsigned int Dimension = 3;
	typedef float CoordType;
	typedef itk::Image< CoordType, Dimension >         ImageType;
	typedef itk::ImageFileReader<ImageType> ReaderType;
	typedef itk::VariableLengthVector<CoordType> MeasurementVectorType;
	typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;
	typedef fs::Surface< CoordType, Dimension> SurfaceType;
	typedef itk::Statistics::WeightedCentroidKdTreeGenerator< SampleType > TreeGeneratorType;
	typedef TreeGeneratorType::KdTreeType TreeType;
	typedef itk::Statistics::KdTreeBasedKmeansEstimator<TreeType> EstimatorType;
	typedef itk::Statistics::DistanceToCentroidMembershipFunction< MeasurementVectorType >	MembershipFunctionType;
	typedef MembershipFunctionType::Pointer                      MembershipFunctionPointer;
	typedef itk::Statistics::MinimumDecisionRule DecisionRuleType;
	typedef itk::Statistics::SampleClassifierFilter< SampleType > ClassifierType;
	typedef ClassifierType::ClassLabelVectorObjectType               ClassLabelVectorObjectType;
	typedef ClassifierType::ClassLabelVectorType                     ClassLabelVectorType;
	typedef ClassifierType::MembershipFunctionVectorObjectType       MembershipFunctionVectorObjectType;
	typedef ClassifierType::MembershipFunctionVectorType             MembershipFunctionVectorType;
	typedef SurfaceType::TriangleType 				TriangleType;
  	
	GetPot cl(narg, const_cast<char**>(arg));
	if(cl.size()==1 || cl.search(2,"--help","-h"))
	{
		std::cout<<"Usage: " << std::endl;
		std::cout<< arg[0] << " -s surface -i image  -c numClusters -d deep -o outputImage"  << std::endl;   
		return -1;
	}
	const char *surfFilename= cl.follow ("", "-s");
	const char *imageFilename = cl.follow ("", "-i");
	const char *outputImageFilename = cl.follow ("", "-o");
	int numClusters = cl.follow(10,"-c");
	int deep= cl.follow(10,"-d");
	
	MRI_SURFACE *surf;
	surf = MRISread(surfFilename);

	SurfaceType::Pointer surface =  SurfaceType::New();
	surface->Load(&*surf);

 	MRI *imageFS =  MRIread(imageFilename) ;

	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(imageFilename);
	reader->Update();
	ImageType::Pointer image = reader->GetOutput();

	SampleType::Pointer sample = SampleType::New();
	sample->SetMeasurementVectorSize( deep*2);

	std::vector<SurfaceType::PointType>  points;
	ImageType::Pointer output = ImageType::New();
	output->SetRegions(image->GetLargestPossibleRegion());
	output->Allocate();
	std::cout << image->GetDirection() << std::endl;
	ImageType::DirectionType direction =  image->GetDirection();
	/*for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			direction[i][j]=-direction[i][j];
		}
	}
*/
	//image->SetDirection(image->GetInverseDirection());
  	//std::cout << image->GetInverseDirection()  << std::endl;
	//image->SetDirection(direction);
	
	output->SetDirection(direction);

	output->SetOrigin(image->GetOrigin());
  	output->SetSpacing(image->GetSpacing());
	
for(SurfaceType::CellsContainerConstIterator itCells = surface->GetCells()->Begin();itCells != surface->GetCells()->End();itCells++)
	{
		SurfaceType::CellType::PointIdIterator pointsIt = itCells.Value()->PointIdsBegin();
		SurfaceType::PointType p1 = surface->GetPoint(*pointsIt);
		SurfaceType::PointType p2 = surface->GetPoint(*(++pointsIt));
		SurfaceType::PointType p3 = surface->GetPoint(*(++pointsIt));
		SurfaceType::PointType edge1 = p1-p2;
		SurfaceType::PointType edge2 = p1-p3;
		std::cout << p1 << std::endl;	
		vnl_vector<CoordType> normal = vnl_cross_3d(edge1.GetVnlVector(), edge2.GetVnlVector());
		normal=normal.normalize();
		SurfaceType::PointType mean ;
		mean.Fill(0.0);
		for(int i=0;i<3;i++)
			mean[i]= (p1[i]+p2[i]+p3[i])/3.0;
	
		ImageType::IndexType index; 
	
		MeasurementVectorType mv;
		mv.SetSize(deep*2);
	
		for( int d=-deep;d<deep;d++)
		{
			SurfaceType::PointType pt;
			for(int i=0;i<3;i++)
				pt[i]= mean[i]+d*normal[i];

			/*image->TransformPhysicalPointToIndex(pt, index);
			std::cout << index << std::endl;

			output->TransformPhysicalPointToIndex(pt, index);
			std::cout << index << std::endl;
			*/
			double x,y,z;
			MRISsurfaceRASToVoxel(surf, imageFS, pt[0], pt[1],pt[2], &x,&y,&z);
			//if (output->TransformPhysicalPointToIndex(pt, index))
			index[0]=x;
			index[1]=y;
			index[2]=z;
			mv[d+deep]= image->GetPixel(index);
		}
		sample->PushBack( mv );
		points.push_back(mean);
	}

	//output->SetDirection(direction);
  	TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();
	treeGenerator->SetSample( sample );
	treeGenerator->SetBucketSize( 160 );
	treeGenerator->Update();

	EstimatorType::Pointer estimator = EstimatorType::New();
	EstimatorType::ParametersType initialMeans(deep*numClusters);
	for(unsigned int c=0;c<numClusters;c++)
	{
		for(unsigned int i=0;i<deep;i++)
			initialMeans[c*deep+i] = sample->GetMeasurementVector(c*(sample->Size()/numClusters))[i]; 
		
	}
	estimator->SetParameters( initialMeans );
	estimator->SetKdTree( treeGenerator->GetOutput() );
	estimator->SetMaximumIteration( 200 );
	estimator->SetCentroidPositionChangesThreshold(0.0);
	estimator->StartOptimization();

	EstimatorType::ParametersType estimatedMeans = estimator->GetParameters();

	for ( unsigned int i = 0 ; i < numClusters ; ++i )
	{
		std::cout << "cluster[" << i << "] " << std::endl;
		std::cout << "    estimated mean : " << estimatedMeans[i] << std::endl;
	}

	DecisionRuleType::Pointer decisionRule = DecisionRuleType::New();
	ClassifierType::Pointer classifier = ClassifierType::New();
	classifier->SetDecisionRule(decisionRule);
	classifier->SetInput( sample );
	classifier->SetNumberOfClasses( numClusters );

	ClassLabelVectorObjectType::Pointer  classLabelsObject = ClassLabelVectorObjectType::New();
	classifier->SetClassLabels( classLabelsObject );

	ClassLabelVectorType &  classLabelsVector = classLabelsObject->Get();
	for(unsigned int i=0; i<numClusters;i++)
		classLabelsVector.push_back( i );


	MembershipFunctionVectorObjectType::Pointer membershipFunctionsObject =
	MembershipFunctionVectorObjectType::New();
	classifier->SetMembershipFunctions( membershipFunctionsObject );

	MembershipFunctionVectorType &  membershipFunctionsVector = membershipFunctionsObject->Get();

	MembershipFunctionType::CentroidType origin( sample->GetMeasurementVectorSize() );
	int index = 0;
	std::vector<MembershipFunctionPointer> functions;
	std::vector<MeasurementVectorType> clusterMeans;
	for ( unsigned int i = 0 ; i < numClusters ; i++ )
	{
		MembershipFunctionPointer membershipFunction = MembershipFunctionType::New();
		for ( unsigned int j = 0 ; j < sample->GetMeasurementVectorSize(); j++ )
		{
			origin[j] = estimatedMeans[index++];
		}
		membershipFunction->SetCentroid( origin );
		functions.push_back(membershipFunction);
		membershipFunctionsVector.push_back( membershipFunction.GetPointer() );
		MeasurementVectorType mean;
		mean.SetSize(deep*2);
		mean.Fill(0);
		clusterMeans.push_back(mean);
	}
	classifier->Update();

	const ClassifierType::MembershipSampleType* membershipSample = classifier->GetOutput();
	ClassifierType::MembershipSampleType::ConstIterator iter = membershipSample->Begin();
	std::vector<int> norm(numClusters, 0);
	while ( iter != membershipSample->End() )
	{
		clusterMeans[iter.GetClassLabel()]+=iter.GetMeasurementVector();
		norm[iter.GetClassLabel()]++;
		++iter;
	}	
	for (int i=0;i<numClusters;i++)
	{
		std::cout << clusterMeans[i]/norm[i]<< " " << norm[i] << std::endl;
		if(norm[i]>0)
			clusterMeans[i]/=norm[i];
	}	

	iter = membershipSample->Begin();
	int i=0;
	while ( iter != membershipSample->End() )
	{

		//MeasurementVectorType meanIntensityDecay = clusterMeans[iter.GetClassLabel()]; ///norm[iter.GetClassLabel()];
		ImageType::IndexType index;
		//output->TransformPhysicalPointToIndex(points[i], index);

		double x,y,z;
		MRISsurfaceRASToVoxel(surf, imageFS, points[i][0], points[i][1],points[i][2], &x,&y,&z);
		index[0]=x;
		index[1]=y;
		index[2]=z;
		output->SetPixel(index, iter.GetClassLabel());
		
		++iter;
		i++;
	}
	typedef  itk::ImageFileWriter< ImageType  > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(outputImageFilename);
	writer->SetInput(output);
	writer->Update();
}
