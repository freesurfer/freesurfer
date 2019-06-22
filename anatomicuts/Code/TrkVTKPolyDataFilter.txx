#ifndef _TrkVTKPolyDataFilter_txx_
#define _TrkVTKPolyDataFilter_txx_

#include "TrkVTKPolyDataFilter.h"
//#include "vial.h"	// Needs to be included first because of CVS libs
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkFieldData.h>


#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkFieldData.h>
#include <vtkSmartPointer.h>
#include "itkPoint.h"
#include "itkImage.h"
#include "itkContinuousIndex.h"
#include <vnl/vnl_inverse.h>
	template<class TImage>
TrkVTKPolyDataFilter<TImage>::TrkVTKPolyDataFilter()
{
	m_vtk = vtkPolyData::New();
	//this->SetNumberOfInputs (1);
	//unsigned char color[3] = {219,112,147};
	this->m_color =  nullptr;
	this->m_refHeader = 0;
}

template<class TImage>
TrkVTKPolyDataFilter<TImage>::~TrkVTKPolyDataFilter()
{
	//m_vtk->Delete();
}


	template<class TImage>
void TrkVTKPolyDataFilter<TImage>::GenerateData()
{
}

	template<class TImage>
void TrkVTKPolyDataFilter<TImage>::TrkToVTK()
{

	m_vtk->Initialize();
	m_vtk->Allocate();

	vtkUnsignedCharArray* allColors = vtkUnsignedCharArray::New();
	allColors->SetNumberOfComponents (3);

	vtkSmartPointer<vtkPoints> points = vtkPoints::New();
	
	CTrackReader trkreader;
	TRACK_HEADER trkheadin;
	trkreader.Open(this->m_trkFileName.c_str(), &trkheadin);
	int npts= 0;
	int totalPoints = 0;

	//std::vector<int> orientation;
	//orientation.push_back((trkheadin.voxel_order_original[0]  == 'L')?1:-1);
	//orientation.push_back((trkheadin.voxel_order_original[1]  == 'P')?1:-1);
	//orientation.push_back((trkheadin.voxel_order_original[2]  == 'S')?1:-1);
	//std::cout << orientation[0]<< trkheadin.voxel_order_original[0] << std::endl;
	//std::cout << orientation[1]<< trkheadin.voxel_order_original[1] << std::endl;
	//std::cout << orientation[2]<< trkheadin.voxel_order_original[2] << std::endl;
	while (trkreader.GetNextPointCount(&npts))
	{
		float *iraw, *rawpts = new float[npts*3];
		std::vector<float> newpts;

		// Read a streamline from input file
		trkreader.GetNextTrackData(npts, rawpts);

		iraw = rawpts;
		itk::Point<float> pt;
		itk::ContinuousIndex<float,3> index;
		itk::ContinuousIndex<float,4> index4;
		vtkIdType *ids = new vtkIdType [npts];

		for (int ipt=0 ; ipt <npts; ipt++) {
			// Divide by input voxel size and make 0-based to get voxel coords
			for (int k = 0; k < 3; k++) 
			{
				pt[k] = *iraw  ;// / trkheadin.voxel_size[k]  ;//*trkheadin.voxel_size[k];
				index[k]= pt[k]/trkheadin.voxel_size[k]; //*orientation[k];
				index4[k]= index[k];
				iraw++;
			}
			index4[3]=1;
			if (!this->m_refImage.IsNull())
			//if(trkheadin.vox_to_ras[3][3]==0) //not recorded
			{
				//std::cout << " reference iamge"  << std::endl;
				this->m_refImage->TransformContinuousIndexToPhysicalPoint(index,pt);	
			}
			else
			{
				//pt.Fill(0.0);
				//this->m_refImage->TransformContinuousIndexToPhysicalPoint(index,pt);	
//				std::cout << "pt1 "<< pt << std::endl;
				pt.Fill(0.0);
				for (int k1 = 0; k1 < 3; k1++) 
				{
					for (int k2 = 0; k2 < 4; k2++) 
					{
						pt[k1] += index4[k2]*trkheadin.vox_to_ras[k1][k2];

					}
				}
//				std::cout <<"p2 "<< pt << std::endl;

			}

			
			points->InsertPoint (totalPoints, pt[0], pt[1], pt[2] );
			ids[ipt] = totalPoints;
			totalPoints++;
		}

		m_vtk->InsertNextCell (VTK_POLY_LINE, npts, ids);
		if ( this->m_color != nullptr)
		{
			allColors->InsertNextTuple3 ( this->m_color[0], this->m_color[1], this->m_color[2] );
		}
		delete [] ids;
	}

	m_vtk->SetPoints ( points );
	if(allColors->GetSize()>0)
	{	
		m_vtk->GetCellData()->SetScalars ( allColors );
	}
	points->Delete();

}

	template<class TImage>
void TrkVTKPolyDataFilter<TImage>::VTKToTrk(std::string outputName)
{

	CTrackWriter trkwriter;
	TRACK_HEADER trkheadout;
	
	if (this->m_refHeader ==0)
	{
		trkheadout.Initialize();

		//std::cout << "hola " << std::endl;
		for (int i=0; i<3 ; i++)
		{
			trkheadout.origin[i] =m_refImage->GetOrigin()[i]; 
			trkheadout.voxel_size[i] = m_refImage->GetSpacing()[i];
			trkheadout.dim[i] = m_refImage->GetLargestPossibleRegion().GetSize()[i];

			for(int j=0;j<3;j++)
			{
				if(i==j)
					trkheadout.vox_to_ras[i][j]= m_refImage->GetDirection()[i][j]*m_refImage->GetSpacing()[i];
				else
					trkheadout.vox_to_ras[i][j]= m_refImage->GetDirection()[i][j];
			}
			trkheadout.vox_to_ras[i][3]= m_refImage->GetOrigin()[i];//*m_refImage->GetSpacing()[i];
			
		}
		//std::cout << m_refImage->GetLargestPossibleRegion().GetIndex()[0]<<std::endl;
		//std::cout << m_refImage->GetOrigin()[0]<< std::endl;
		trkheadout.vox_to_ras[3][3]=1;
		//std::cout << std::endl;
		trkheadout.image_orientation_patient[0] =   trkheadout.vox_to_ras[0][0] / trkheadout.voxel_size[0];
		trkheadout.image_orientation_patient[1] =   trkheadout.vox_to_ras[1][0] / trkheadout.voxel_size[0];
		trkheadout.image_orientation_patient[2] =  trkheadout.vox_to_ras[2][0] / trkheadout.voxel_size[0];
		trkheadout.image_orientation_patient[3] =   trkheadout.vox_to_ras[0][1] / trkheadout.voxel_size[1];
		trkheadout.image_orientation_patient[4] =   trkheadout.vox_to_ras[1][1] / trkheadout.voxel_size[1];
		trkheadout.image_orientation_patient[5] =  trkheadout.vox_to_ras[2][1] / trkheadout.voxel_size[1];
	}
	else
	{
		//trkwriter.Initialize(outputName.c_str(), m_refHeader);
		trkheadout = *this->m_refHeader;
	}
	if( this->m_color != nullptr ) 
	{
		trkheadout.reserved[0] = 'C';
		for (int i=1;i<4;i++)
		{
			trkheadout.reserved[i] = this->m_color[i-1] ;
			//std::c	out << this->m_color[i-1] << std::endl;
		}
	}
	/*float hola[4][4]; //= new float[4][4]();

	for (int i=0; i<4 ; i++)
	{	
		for(int j=0;j<4;j++)
		{
			hola[i][j]=i+j; // m_refImage->GetDirection()[i][j]*m_refImage->GetSpacing()[i];
			std::cout << trkheadout.vox_to_ras[i][j]<< " ";
		}
	}
	trkheadout.vox_to_ras=hola;
*/

	trkwriter.Initialize(outputName.c_str(), trkheadout);
	trkwriter.UpdateHeader(trkheadout);
	vtkCellArray *lines = m_vtk->GetLines();
	lines->InitTraversal();
	vtkIdType pointCount, *pointBuf;
	vnl_matrix<float> vox_to_ras = vnl_matrix<float>(4,4);
	for (int k1 = 0; k1 < 4; k1++) 
		for (int k2 = 0; k2 < 4; k2++) 
			vox_to_ras(k1,k2)= trkheadout.vox_to_ras[k1][k2];


	vnl_matrix<float> ras_to_vox = 	vnl_inverse(vox_to_ras);
	//std::cout << vox_to_ras<< std::endl;
	//std::cout << ras_to_vox << std::endl;
	while ( lines->GetNextCell(pointCount, pointBuf) )
	{

		std::vector<float> points( pointCount*3,0);
		int n=0;
		for (vtkIdType k=0; k<pointCount; k++)
		{
			double *pt = m_vtk->GetPoint ( pointBuf[k] );

			itk::Point<float> pt2;
			itk::ContinuousIndex<float,3> index;
			itk::ContinuousIndex<float,4> point4;
			for (int k = 0; k < 3; k++) {
				pt2[k]= pt[k];
				point4[k]=pt[k];
			}
			point4[3]=1;
			
			if(! this->m_refImage.IsNull())
			//if(this->m_refHeader==0) //not recorded
			{
				this->m_refImage->TransformPhysicalPointToContinuousIndex(pt2,index);	
				for (unsigned int i=0; i<3; i++)
				{
					points [n] = index[i]* trkheadout.voxel_size[i];
					n++;
				}
			}
			else
			{
				pt2.Fill(0.0);	
				for (int k1 = 0; k1 < 3; k1++) 
				{
					for (int k2 = 0; k2 < 4; k2++) 
					{
						pt2[k1] += point4[k2]*ras_to_vox[k1][k2];
					}
				}


				for (unsigned int i=0; i<3; i++)
				{
					points [n] = pt2[i]*trkheadout.voxel_size[i];
					n++;
				}
			}
		}
		trkwriter.WriteNextTrack(pointCount, &(points.at(0)));
	}


	trkwriter.Close();
}

/*
   vtkSmartPointer<vtkIntArray> intArrayCellData  = vtkIntArray::New();
   typedef typename MeshType::CellDataContainer CellDataContainer;
   if(  input->GetCellData() != 0 )
   {

   typename CellDataContainer::ConstIterator cellData = input->GetCellData()->Begin();
   int i=0;
   for(;cellData!= input->GetCellData()->End(); ++cellData)
   {
   intArrayCellData->InsertValue(i,cellData.Value());
   i++;
   }

   vtkSmartPointer<vtkFieldData> fieldData = vtkFieldData::New();
   fieldData->AddArray(intArrayCellData);
   m_Output->SetFieldData(fieldData);

   }*/
/*


#endif
#ifndef _TrkVTKPolyDataFilter_txx_
#define _TrkVTKPolyDataFilter_txx_

#include "TrkVTKPolyDataFilter.h"
//#include "vial.h"	// Needs to be included first because of CVS libs
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkFieldData.h>


#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkFieldData.h>
#include <vtkSmartPointer.h>
#include "itkPoint.h"
#include "itkImage.h"
#include "itkContinuousIndex.h"

	template<class TImage>
TrkVTKPolyDataFilter<TImage>::TrkVTKPolyDataFilter()
{
	m_vtk = vtkPolyData::New();
	this->SetNumberOfInputs (1);
	unsigned char color[3] = {219,112,147};
	this->m_color =  color;
}

template<class TImage>
TrkVTKPolyDataFilter<TImage>::~TrkVTKPolyDataFilter()
{
	//m_vtk->Delete();
}


	template<class TImage>
void TrkVTKPolyDataFilter<TImage>::GenerateData()
{
}

	template<class TImage>
void TrkVTKPolyDataFilter<TImage>::TrkToVTK()
{

	m_vtk->Initialize();
	m_vtk->Allocate();

	vtkUnsignedCharArray* allColors = vtkUnsignedCharArray::New();
	allColors->SetNumberOfComponents (3);

	vtkSmartPointer<vtkPoints> points = vtkPoints::New();
	
	CTrackReader trkreader;
	TRACK_HEADER trkheadin;
	trkreader.Open(this->m_trkFileName.c_str(), &trkheadin);
	int npts= 0;
	int totalPoints = 0;

	std::vector<int> orientation;
	orientation.push_back((trkheadin.voxel_order_original[0]  == 'L')?1:-1);
	orientation.push_back((trkheadin.voxel_order_original[1]  == 'P')?1:-1);
	orientation.push_back((trkheadin.voxel_order_original[2]  == 'S')?1:-1);
	//std::cout << orientation[0]<< trkheadin.voxel_order_original[0] << std::endl;
	//std::cout << orientation[1]<< trkheadin.voxel_order_original[1] << std::endl;
	//std::cout << orientation[2]<< trkheadin.voxel_order_original[2] << std::endl;
	while (trkreader.GetNextPointCount(&npts))
	{
		float *iraw, *rawpts = new float[npts*3];
		std::vector<float> newpts;

		// Read a streamline from input file
		trkreader.GetNextTrackData(npts, rawpts);

		iraw = rawpts;
		itk::Point<float> pt;
		itk::ContinuousIndex<float,3> index;
		vtkIdType *ids = new vtkIdType [npts];

		for (int ipt=0 ; ipt <npts; ipt++) {
			// Divide by input voxel size and make 0-based to get voxel coords
		

			for (int k = 0; k < 3; k++) {
				pt[k] = *iraw  ;// / trkheadin.voxel_size[k]  ;// *trkheadin.voxel_size[k];
				index[k]= pt[k]/trkheadin.voxel_size[k];
				iraw++;
			}
			if (!this->m_refImage.IsNull())
			{
				this->m_refImage->TransformContinuousIndexToPhysicalPoint(index,pt);	
			}
			points->InsertPoint (totalPoints, pt[0], pt[1], pt[2] );
			ids[ipt] = totalPoints;
			totalPoints++;
		}

		m_vtk->InsertNextCell (VTK_POLY_LINE, npts, ids);
		if ( this->m_color != NULL)
		{
			allColors->InsertNextTuple3 ( this->m_color[0], this->m_color[1], this->m_color[2] );
		}
		delete [] ids;
	}

	m_vtk->SetPoints ( points );
	m_vtk->GetCellData()->SetScalars ( allColors );
	points->Delete();

}

	template<class TImage>
void TrkVTKPolyDataFilter<TImage>::VTKToTrk(std::string outputName)
{

	CTrackWriter trkwriter;
	TRACK_HEADER trkheadout;
	
	if (! m_refImage.IsNull())
	{
		std::cout << "hola " << std::endl;
		for (int i=0; i<3 ; i++)
		{
			trkheadout.origin[i] =m_refImage->GetOrigin()[i]; 
			trkheadout.voxel_size[i] = m_refImage->GetSpacing()[i];
			trkheadout.dim[i] = m_refImage->GetLargestPossibleRegion().GetSize()[i];

			for(int j=0;j<3;j++)
			{
				trkheadout.vox_to_ras[i][j]= m_refImage->GetDirection()[i][j];
			}
		}
		trkheadout.image_orientation_patient[0] =  - trkheadout.vox_to_ras[0][0] / trkheadout.voxel_size[0];
		trkheadout.image_orientation_patient[1] =  - trkheadout.vox_to_ras[1][0] / trkheadout.voxel_size[0];
		trkheadout.image_orientation_patient[2] =  trkheadout.vox_to_ras[2][0] / trkheadout.voxel_size[0];
		trkheadout.image_orientation_patient[3] =  - trkheadout.vox_to_ras[0][1] / trkheadout.voxel_size[1];
		trkheadout.image_orientation_patient[4] =  - trkheadout.vox_to_ras[1][1] / trkheadout.voxel_size[1];
		trkheadout.image_orientation_patient[5] =  trkheadout.vox_to_ras[2][1] / trkheadout.voxel_size[1];

		trkheadout.reserved[0] = 'C';
		for (int i=1;i<4;i++)
			trkheadout.reserved[i] = this->m_color[i-1] ;
		
		trkwriter.Initialize(outputName.c_str(), trkheadout);
		trkwriter.UpdateHeader(trkheadout);
	}
	else
	{
		std::cout << "pucha" << std::endl;
		trkwriter.Initialize(outputName.c_str(), m_refHeader);
	}
	vtkCellArray *lines = m_vtk->GetLines();
	lines->InitTraversal();
	vtkIdType pointCount, *pointBuf;
	while ( lines->GetNextCell(pointCount, pointBuf) )
	{

		std::vector<float> points( pointCount*3,0);
		int n=0;
		for (vtkIdType k=0; k<pointCount; k++)
		{
			double *pt = m_vtk->GetPoint ( pointBuf[k] );

			itk::Point<float> pt2;
			itk::ContinuousIndex<float,3> index;
			for (int k = 0; k < 3; k++) {
				pt2[k]= pt[k];
			}
			if(! this->m_refImage.IsNull())
			{
				this->m_refImage->TransformPhysicalPointToContinuousIndex(pt2,index);	
				for (unsigned int i=0; i<3; i++)
				{
					points [n] = index[i]* trkheadout.voxel_size[i];
					n++;
				}
			}
			else
			{
				for (unsigned int i=0; i<3; i++)
				{
					points [n] = pt[i];
					n++;
				}
			}
		}
		trkwriter.WriteNextTrack(pointCount, &(points.at(0)));
	}


	trkwriter.Close();
}

*/

#endif
