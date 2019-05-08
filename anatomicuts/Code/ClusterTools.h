#ifndef _ClusterTools_h__
#define _ClusterTools_h_

#include "itkMeshToMeshFilter.h"
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include "vtkSplineFilter.h"
#include "OrientationPlanesFromParcellationFilter.h"
#include "LabelPerPointVariableLengthVector.h"
#include "LabelsEntropyAndIntersectionMembershipFunction.h"

#include "itkDefaultStaticMeshTraits.h"

#include "EuclideanMembershipFunction.h"
#include "LabelPerPointVariableLengthVector.h"
#include "colortab.h"
#include "fsenv.h"

#include "TrkVTKPolyDataFilter.txx"

#include "itkPolylineCell.h"
#include "vtkCellData.h"

#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkFieldData.h>
#include "PolylineMeshToVTKPolyDataFilter.h"
#include "VTKPolyDataToPolylineMeshFilter.h"

using namespace itk;
enum DirectionsType{ALL=0, DIAGONAL=1, STRAIGHT=2};

template <class TColorMesh, class TImage, class THistogramMesh>
class ClusterTools:	public LightObject //<TColorMesh, TColorMesh>
{
	public:
		typedef ClusterTools Self;
		typedef MeshToMeshFilter<TColorMesh, TColorMesh>       Superclass;
		typedef SmartPointer<Self>                              Pointer;
		typedef SmartPointer<const Self>                        ConstPointer;

		itkNewMacro  (Self);
		itkTypeMacro (ClusterTools, LightObject); // MeshToMeshFilter);

		typedef TColorMesh                         ColorMeshType;
		typedef typename ColorMeshType::Pointer  ColorMeshPointer;
		typedef typename ColorMeshType::PointType  ColorPointType;
		typedef typename ColorMeshType::PixelType  ColorPixelType;


		typedef typename  ColorMeshType::CellType              CellType;
		typedef typename CellType::CellAutoPointer       CellAutoPointer;
		typedef PolylineMeshToVTKPolyDataFilter<ColorMeshType> VTKConverterType;
		typedef VTKPolyDataToPolylineMeshFilter<ColorMeshType> MeshConverterType;
		
		typedef THistogramMesh                         HistogramMeshType;
		typedef typename HistogramMeshType::Pointer  HistogramMeshPointer;
		typedef typename HistogramMeshType::PointType  HistogramPointType;
		typedef typename HistogramMeshType::PixelType  HistogramPixelType;
		typedef typename HistogramMeshType::CellPixelType HistogramDataType;
		typedef LabelPerPointVariableLengthVector<ColorPixelType , HistogramMeshType> MeasurementVectorType;	
		typedef LabelsEntropyAndIntersectionMembershipFunction<MeasurementVectorType>  MembershipFunctionType;	
		//typedef EuclideanMembershipFunction<MeasurementVectorType> MembershipFunctionType;

		typedef typename ColorMeshType::PointsContainerPointer PointsContainerPointer;

		typedef typename ColorMeshType::PointsContainer 		PointsContainer;

		typedef TImage ImageType;
		typedef typename ImageType::Pointer ImagePointer;
		typedef typename ImageType::IndexType IndexType;

		void GetPolyDatas(std::vector<std::string> files, std::vector<vtkSmartPointer<vtkPolyData>>* polydatas, ImagePointer image);
		std::vector<ColorMeshPointer>* FixSampleClusters(std::vector<vtkSmartPointer<vtkPolyData>> p , int i);
		std::vector<ColorMeshPointer>* PolydataToMesh(std::vector<vtkSmartPointer<vtkPolyData>> p );
		std::vector<HistogramMeshPointer>* ColorMeshToHistogramMesh(std::vector<ColorMeshPointer> basicMeshes, ImagePointer segmentation, bool removeInterHemispheric);
		void  SetDirectionalNeighbors(std::vector<HistogramMeshPointer>* meshes, ImagePointer segmentation, std::vector<itk::Vector<float>> direcciones, bool symmetry);
	
		int GetAverageStreamline(ColorMeshPointer mesh);
		float  GetStandardDeviation(HistogramMeshPointer mesh, int averageStreamlineIndex);
		float GetDistance(HistogramMeshPointer mesh, int index_i, int index_j); 
		
			//typedef std::vector<int>                  PointDataType;
		
		
		void SaveMesh(ColorMeshPointer mesh, ImagePointer image,std::string outputFilename, std::string refFiber);


		int SymmetricLabelId(int);
		std::vector<itk::Vector<float>>* GetDirections(DirectionsType dir);
		
	protected:
		ClusterTools(){}
		~ClusterTools() {}


	private:
		ClusterTools (const Self&);
		void operator=(const Self&);
};


#include "ClusterTools.txx"
#endif
