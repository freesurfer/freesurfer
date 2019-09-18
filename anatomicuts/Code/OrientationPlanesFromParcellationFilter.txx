#ifndef __OrientationPlanesFromParcellationFilter_txx
#define __OrientationPlanesFromParcellationFilter_txx

#include "OrientationPlanesFromParcellationFilter.h"

template<class TInputImage,class TOutputImage>
void OrientationPlanesFromParcellationFilter<TInputImage,TOutputImage>::GenerateData()
{

	vtkSmartPointer<vtkPoints> pointsCC_3V = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPoints> pointsCC1 = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPoints> pointsCC5 = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPoints> rhThalamus= vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPoints> lhThalamus = vtkSmartPointer<vtkPoints>::New();

	typename InputImageType::SizeType radius;
	radius.Fill(1); 
	itk::ConstNeighborhoodIterator<InputImageType> iter(radius, this->GetInput(),this->GetInput()->GetLargestPossibleRegion());

	for(;!iter.IsAtEnd();++iter)
	{
		typename InputImageType::PixelType pixel = iter.GetCenterPixel();
		typename InputImageType::PointType point;
		this->GetInput()->TransformIndexToPhysicalPoint(iter.GetIndex(),point);

		if(pixel == 14 || pixel == 254 || pixel ==253 || pixel ==252 ||   (pixel ==3026 && m_baby)|| (pixel ==4026 && m_baby) || (pixel== 3010 && m_baby)|| (pixel== 4010 && m_baby) || (pixel ==175 && m_baby)) // || pixel == 255 || pixel == 251)
		{
			pointsCC_3V->InsertNextPoint(point[0],point[1], point[2]);
		}
		if ( (pixel == 251 && !m_baby) || (pixel ==3026 && m_baby)|| (pixel ==4026 && m_baby)  )
		{
			pointsCC1->InsertNextPoint(point[0],point[1], point[2]);
		}
		if ((pixel==255 && !m_baby ) || (pixel== 3010 && m_baby)|| (pixel== 4010 && m_baby)  )
		{
			pointsCC5->InsertNextPoint(point[0],point[1], point[2]);
		}
		if(pixel==49 || pixel ==48)
		{
			rhThalamus->InsertNextPoint(point[0],point[1], point[2]);
		}
		if(pixel==10 || pixel ==9)
		{
			lhThalamus->InsertNextPoint(point[0],point[1], point[2]);
		}

	}
/*	std::cout <<  pointsCC_3V->GetNumberOfPoints() << std::endl;
	std::cout <<  pointsCC1->GetNumberOfPoints() << std::endl;
	std::cout <<  pointsCC5->GetNumberOfPoints() << std::endl;
	std::cout <<  rhThalamus->GetNumberOfPoints() << std::endl;
	std::cout <<  lhThalamus->GetNumberOfPoints() << std::endl;
*/
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(pointsCC_3V);

	// Compute the center of mass
	vtkSmartPointer<vtkCenterOfMass> centerOfMassFilter =	vtkSmartPointer<vtkCenterOfMass>::New();
#if VTK_MAJOR_VERSION > 5
	centerOfMassFilter->SetInputData(polydata);
#else
	centerOfMassFilter->SetInput(polydata);
#endif
	centerOfMassFilter->SetUseScalarsAsWeights(false);
	centerOfMassFilter->Update();
	
	double center[3];
	centerOfMassFilter->GetCenter(center);

	double **a = new double *[3];     
	double *eigval = new double [3];
	double **eigvec = new double *[3];

	// Clear the upper triangular region (and btw, allocate the eigenvecs as well)
	for (unsigned int i = 0; i <3; i++)
	{
		eigvec[i] = new double[3];
		a[i] = new double[3];
		eigval[i]=0.0;
		for (unsigned int j = 0; j < 3; j++)
		{
			a[i][j] = 0.0;
			eigvec[i][j]=0.0;
		}
	}

	double weightSum = 0;
	for(vtkIdType pointId = 0; pointId < pointsCC_3V->GetNumberOfPoints(); pointId++ )
	{
		double x[3];
		double xp[3];
		pointsCC_3V->GetPoint(pointId, x);

		for(unsigned int j = 0; j < 3; j++ )
		{
			xp[j] = x[j] - center[j];

			for (unsigned int i = j; i < 3; i++)
			{
				double value = xp[j] * xp[i];
				a[j][i] += value;
				a[i][j] += value;
			}
			weightSum += 1;
		}
	}

	// Divide by N-1 for an unbiased estimate
	for(unsigned int  j = 0; j < 3; j++ )
	{
		for(unsigned int i = 0; i < 3; i++)
		{
			a[j][i] /= weightSum;
		}
	}

	// Extract eigenvectors from covariance matrix
	vtkMath::Jacobi(a,eigval,eigvec);

	double normal[3];
	normal[0] = eigvec[0][2];
	normal[1] = eigvec[1][2];
	normal[2] = eigvec[2][2];
	std::cout << normal[0] << normal[1] << normal[2] <<std::endl;

	//orient the normal of the sagital plane from left to right
	centerOfMassFilter =	vtkSmartPointer<vtkCenterOfMass>::New();
	polydata->SetPoints(lhThalamus);
#if VTK_MAJOR_VERSION > 5
	centerOfMassFilter->SetInputData(polydata);
#else
	centerOfMassFilter->SetInput(polydata);
#endif
	centerOfMassFilter->SetUseScalarsAsWeights(false);
	centerOfMassFilter->Update();
	
	double lhThalamusCenter[3];
	centerOfMassFilter->GetCenter(lhThalamusCenter);

	centerOfMassFilter =	vtkSmartPointer<vtkCenterOfMass>::New();
	polydata->SetPoints(rhThalamus);
#if VTK_MAJOR_VERSION > 5
	centerOfMassFilter->SetInputData(polydata);
#else
	centerOfMassFilter->SetInput(polydata);
#endif
	centerOfMassFilter->SetUseScalarsAsWeights(false);
	centerOfMassFilter->Update();
	
	double rhThalamusCenter[3];
	centerOfMassFilter->GetCenter(rhThalamusCenter);
	double directionLhRh[3];
	int innerProduct=0;
	for(int i=0;i<3;i++)
	{
		directionLhRh[i]= rhThalamusCenter[i]- lhThalamusCenter[i];
		innerProduct+=directionLhRh[i]*normal[i];
	}
	if(innerProduct <0)
	{
		for(int i=0;i<3;i++)
		{
			directionLhRh[i]= -directionLhRh[i];		
		}
	}

	//find axial plane
	centerOfMassFilter =	vtkSmartPointer<vtkCenterOfMass>::New();
	polydata->SetPoints(pointsCC1);
#if VTK_MAJOR_VERSION > 5
	centerOfMassFilter->SetInputData(polydata);
#else
	centerOfMassFilter->SetInput(polydata);
#endif
	centerOfMassFilter->SetUseScalarsAsWeights(false);
	centerOfMassFilter->Update();
	
	double centerCC1[3];
	centerOfMassFilter->GetCenter(centerCC1);

	centerOfMassFilter =	vtkSmartPointer<vtkCenterOfMass>::New();
	polydata->SetPoints(pointsCC5);
#if VTK_MAJOR_VERSION > 5
	centerOfMassFilter->SetInputData(polydata);
#else
	centerOfMassFilter->SetInput(polydata);
#endif
	centerOfMassFilter->SetUseScalarsAsWeights(false);
	centerOfMassFilter->Update();
	
	double centerCC5[3];
	centerOfMassFilter->GetCenter(centerCC5);

	double projCC1[3], projCC5[3];

	vtkPlane::ProjectPoint(centerCC1,center,normal, projCC1);
	vtkPlane::ProjectPoint(centerCC5,center,normal, projCC5);
		
	double axialCenter[3], thirdPoint[3];
	for(unsigned int i=0;i<3;i++)
	{
		axialCenter[i] = (projCC1[i]+projCC5[i])/2.0;
		thirdPoint[i] = axialCenter[i]+ normal[i]  *10;
	}
	
	vtkSmartPointer<vtkPlaneSource> axial = vtkSmartPointer<vtkPlaneSource>::New();
	axial->SetPoint1(projCC1);	
	axial->SetPoint2(projCC5);	
	//axial->SetPoint3(thirdPoint);
	axial->SetOrigin(thirdPoint);
	axial->Update();

	this->GetOutput()->SetRegions(this->GetInput()->GetLargestPossibleRegion());
	this->GetOutput()->Allocate();
	this->GetOutput()->FillBuffer( 0.0 );

	itk::ConstNeighborhoodIterator<InputImageType> iter2(radius, this->GetInput(),this->GetInput()->GetLargestPossibleRegion());
	for(;!iter2.IsAtEnd();++iter2)
	{
		typename InputImageType::IndexType index = iter2.GetIndex();
		typename InputImageType::PointType point;
		this->GetInput()->TransformIndexToPhysicalPoint(index,point);
		double proj[3];
		double pt[3];
		for(int i=0;i<3;i++)
			pt[i]=point[i];
		
		vtkPlane::ProjectPoint(pt,center,normal, proj);

		for(int i=0;i<3;i++)
			point[i]=proj[i];
		if (this->GetOutput()->TransformPhysicalPointToIndex(point, index))
		{
			this->GetOutput()->SetPixel(index, 255);
		}

		vtkPlane::ProjectPoint(pt,axialCenter,axial->GetNormal(), proj);
		for(int i=0;i<3;i++)
			point[i]=proj[i];
		if (this->GetOutput()->TransformPhysicalPointToIndex(point, index))
		{
			this->GetOutput()->SetPixel(index, 155);
		}
	}
	for(int i=0;i<3;i++)
	{
		m_LeftRight[i]=normal[i];
		m_UpDown=axial->GetNormal()[i];

	}	
	m_FrontBack.SetVnlVector(vnl_cross_3d(m_LeftRight.GetVnlVector(), m_UpDown.GetVnlVector()));
	m_FrontBack.Normalize();
	m_UpDown.Normalize();
	m_LeftRight.Normalize();
	delete[] a;
	delete[] eigval;
	delete[] eigvec;
}
#endif
