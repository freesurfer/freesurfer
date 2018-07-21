#ifndef __TrkVTKPolyDataFilter_h
#define __TrkVTKPolyDataFilter_h
#include "TrackIO.h"

#include "itkProcessObject.h"
#include "itkPolylineCell.h"
#include <vtkSmartPointer.h>

#include "vtkCellData.h"
#include "itkImageFileReader.h"

//class vtkPolyData;

template <class TImage>
class TrkVTKPolyDataFilter : public itk::ProcessObject
{
	public:
		typedef TrkVTKPolyDataFilter Self;
		typedef itk::ProcessObject      Superclass;
		typedef itk::SmartPointer<Self>       Pointer;
		typedef itk::SmartPointer<const Self> ConstPointer;

		itkNewMacro  (Self);
		itkTypeMacro (TrkVTKPolyDataFilter, itk::ProcessObject);

		typedef TImage    ImageType;
		typedef typename ImageType::Pointer ImagePointer;

    		typedef typename itk::ImageFileReader<ImageType> ImageReaderType;

		void SetInput (vtkSmartPointer<vtkPolyData> inputVTK)
		{
			this->m_vtk = inputVTK;
		}

		vtkSmartPointer<vtkPolyData> GetOutputPolyData (void)
		{ return m_vtk; }

		virtual void Update (void)
		{ this->GenerateData(); }
		void SetColor ( unsigned char* color)
		{
			this->m_color = color;
		}
		void SetReferenceTrack(std::string trackref)
		{
			CTrackReader trkreader;
			m_refHeader = new TRACK_HEADER();
			trkreader.Open(trackref.c_str(),this->m_refHeader);
		
	
		}
		void SetReferenceImage(ImagePointer image)
		{
			this->m_refImage = image;
		}
		void SetTrkFileName(std::string filename)
		{
			this->m_trkFileName = filename;
		}
			
		void TrkToVTK();
		void VTKToTrk(std::string outputName);
		TrkVTKPolyDataFilter();
		~TrkVTKPolyDataFilter();
	protected:

		virtual void GenerateData (void);

	private:
		TrkVTKPolyDataFilter (const Self&);
		void operator= (const Self&);
		unsigned char* m_color;
		vtkSmartPointer<vtkPolyData> m_vtk;
		std::string m_trkFileName;
		ImagePointer m_refImage;
		TRACK_HEADER* m_refHeader;
};
#endif
/*
#ifndef ITK_MANUAL_INSTANTIATION
#include "TrkVTKPolyDataFilter.txx"
#endif

#endif
#ifndef _TrkVTKPolyDataFilter_h_
#define _TrkVTKPolyDataFilter_h_
#include "TrackIO.h"

#include "itkProcessObject.h"
#include "itkPolylineCell.h"
#include <vtkSmartPointer.h>

#include "vtkCellData.h"
#include "itkImageFileReader.h"

class vtkPolyData;

template <class TImage>
class TrkVTKPolyDataFilter : public itk::ProcessObject
{
	public:
		typedef TrkVTKPolyDataFilter Self;
		typedef itk::ProcessObject      Superclass;
		typedef itk::SmartPointer<Self>       Pointer;
		typedef itk::SmartPointer<const Self> ConstPointer;

		itkNewMacro  (Self);
		itkTypeMacro (TrkVTKPolyDataFilter, itk::ProcessObject);

		typedef TImage    ImageType;
		typedef typename ImageType::Pointer ImagePointer;

    		typedef typename itk::ImageFileReader<ImageType> ImageReaderType;

		void SetInput (vtkSmartPointer<vtkPolyData> inputVTK)
		{
			this->m_vtk = inputVTK;
		}

		vtkSmartPointer<vtkPolyData> GetOutputPolyData (void)
		{ return m_vtk; }

		virtual void Update (void)
		{ this->GenerateData(); }
		void SetColor ( unsigned char* color)
		{
			this->m_color = color;
		}
		void SetReferenceTrack(std::string trackref)
		{
			CTrackReader trkreader;
			trkreader.Open(trackref.c_str(), &this->m_refHeader);
	
		}
		void SetReferenceImage(ImagePointer image)
		{
			this->m_refImage = image;
		}
		void SetTrkFileName(std::string filename)
		{
			this->m_trkFileName = filename;
		}
			
		void TrkToVTK();
		void VTKToTrk(std::string outputName);
		TrkVTKPolyDataFilter();
		~TrkVTKPolyDataFilter();
	protected:

		virtual void GenerateData (void);

	private:
		TrkVTKPolyDataFilter (const Self&);
		void operator= (const Self&);
		unsigned char* m_color;
		vtkSmartPointer<vtkPolyData> m_vtk;
		std::string m_trkFileName;
		ImagePointer m_refImage;
		TRACK_HEADER m_refHeader;
};


#ifndef ITK_MANUAL_INSTANTIATION
#include "TrkVTKPolyDataFilter.txx"
#endif
*/
