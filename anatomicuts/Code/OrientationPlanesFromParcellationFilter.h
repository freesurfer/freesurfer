#ifndef OrientationPlanesFromParcellationFilter_h
#define OrientationPlanesFromParcellationFilter_h

#include "itkImageToImageFilter.h"
#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "itkConstNeighborhoodIterator.h"
#include "vtkPlane.h"
#include "vtkMath.h"
#include "vtkPolyData.h"
#include "vtkPointSet.h"
#include "vtkCenterOfMass.h"
#include "vtkPlaneSource.h"
#include <vnl/vnl_vector.h>
#include <vnl/vnl_cross.h>

template< class TInputImage , class TOutputImage>
class OrientationPlanesFromParcellationFilter: public itk::ImageToImageFilter<TInputImage,TOutputImage>
{
public:
	using Self= OrientationPlanesFromParcellationFilter;
	using Superclass = itk::ImageToImageFilter<TInputImage, TOutputImage>;

	using Pointer =  itk::SmartPointer<Self>;
	using ConstPointer = itk::SmartPointer  <const Self>;

	using InputImageType = TInputImage;
	using OutputImageType = TOutputImage;
 	
	itkNewMacro(Self);
	using VectorType = itk::Vector<float,3>;
	VectorType GetLeftRight()
	{
		return m_LeftRight;
	}

	VectorType GetUpDown()
	{
		return m_UpDown;
	}

	VectorType GetFrontBack()
	{
		return m_FrontBack;
	}

	void SetBabyMode(bool bb)
	{
		this->m_baby = bb;
	}
	
protected:
	void GenerateData() override;

private:
	bool m_baby=false;
	itk::Vector<float, 3> m_LeftRight;
	itk::Vector<float, 3> m_UpDown;
	itk::Vector<float, 3> m_FrontBack;
};
#include "OrientationPlanesFromParcellationFilter.txx"
#endif

