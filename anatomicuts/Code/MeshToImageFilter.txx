#ifndef _MeshToImageFilter_txx_
#define _MeshToImageFilter_txx_

#include "MeshToImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkGaussianSpatialFunction.h>
#include <algorithm>
#include <iostream>
#include <exception>


using  namespace itk;


template<class TInputMesh, class TOutputImage>
	MeshToImageFilter<TInputMesh,TOutputImage>
::MeshToImageFilter()
{
	this->m_usingLabels =false;
}


// Set the output image size.
template<class TInputMesh, class TOutputImage>
void
	MeshToImageFilter<TInputMesh,TOutputImage>
::SetOutputSize( const SizeType & size )
{
	this->m_OutputRegion.SetSize( size );
}


// Get the output image size.
template<class TInputMesh, class TOutputImage>
const typename MeshToImageFilter<TInputMesh,TOutputImage>
::SizeType &
	MeshToImageFilter<TInputMesh,TOutputImage>
::GetOutputSize()
{
	return this->m_OutputRegion.GetSize();
}


// Set the output image index.
template<class TInputMesh, class TOutputImage>
void
	MeshToImageFilter<TInputMesh,TOutputImage>
::SetOutputIndex( const IndexType & index )
{
	this->m_OutputRegion.SetIndex( index );
}


// Get the output image index.
template<class TInputMesh, class TOutputImage>
const typename MeshToImageFilter<TInputMesh,TOutputImage>
::IndexType &
	MeshToImageFilter<TInputMesh, TOutputImage>
::GetOutputIndex()
{
	return this->m_OutputRegion.GetIndex();
}


// Set the output image spacing.
template <class TInputMesh, class TOutputImage>
void 
	MeshToImageFilter<TInputMesh, TOutputImage>
::SetOutputSpacing( const double* spacing )
{
	SpacingType s( spacing );
	this->SetOutputSpacing( s );
}


// Set the output image origin.
template <class TInputMesh, class TOutputImage>
void 
	MeshToImageFilter<TInputMesh, TOutputImage>
::SetOutputOrigin( const double* origin )
{
	PointType p( origin );
	this->SetOutputOrigin( p );
}

// Helper method to set the output parameters based on this image
template <class TInputMesh, class TOutputImage>
void 
	MeshToImageFilter<TInputMesh, TOutputImage>
::SetOutputParametersFromImage ( const ImageBaseType * image )
{
	if( !image )
	{
		itkExceptionMacro(<< "Cannot use a null image reference");
	}
	this->GetOutput()->CopyInformation(image); 
	this->SetOutputOrigin( image->GetOrigin() );
	this->SetOutputSpacing( image->GetSpacing() );
	this->SetOutputDirection( image->GetDirection() );
	this->SetOutputRegion( image->GetLargestPossibleRegion() );
	SetOutputIndex(image->GetLargestPossibleRegion().GetIndex());
	SetOutputSize(image->GetLargestPossibleRegion().GetSize());
	this->GetOutput()->SetBufferedRegion( this->GetOutput()->GetLargestPossibleRegion() );
	this->GetOutput()->Allocate();
	this->GetOutput()->FillBuffer(0);
}


// Inform pipeline of required output region
template <class TInputMesh, class TOutputImage>
void 
	MeshToImageFilter<TInputMesh,TOutputImage>
::GenerateOutputInformation( void )
{
	// call the superclass' implementation of this method
	Superclass::GenerateOutputInformation();

	// get pointer to the output
	typename OutputImageType::Pointer outputPtr = this->GetOutput();
	if ( !outputPtr )
	{
		return;
	}

	//    outputPtr->SetRegion( m_OutputRegion );    
	outputPtr->SetLargestPossibleRegion( m_OutputRegion );    
	outputPtr->SetRequestedRegion( m_OutputRegion );    
	outputPtr->SetSpacing( m_OutputSpacing );
	outputPtr->SetOrigin( m_OutputOrigin );
	outputPtr->SetDirection( m_OutputDirection );
}



template<class TInputMesh, class TOutputImage>
void
	MeshToImageFilter<TInputMesh,TOutputImage>
::GenerateData()
{
	if( this->m_usingLabels )
		return;
	std::cout << " Generate Data MeshToImagefilter" << std::endl;
	this->GetOutput()->SetBufferedRegion( this->GetOutput()->GetRequestedRegion() );
	this->GetOutput()->Allocate();
	this->GetOutput()->FillBuffer(0);
	//      std::cout <<  " region " << this->GetOutput()->GetRequestedRegion() << std::endl;
	//   this->GetOutput()->FillBuffer(0);
	int outsidePixels = 0;
	IndexType index;
	this->GetOutput()->TransformPhysicalPointToIndex(this->GetInput()->GetPoints()->Begin().Value(),index);
	PointType min=this->GetInput()->GetPoints()->Begin().Value();
	PointType max=this->GetInput()->GetPoints()->Begin().Value();
	//			 max=index;
	for(typename MeshType::PointsContainer::Iterator it=this->GetInput()->GetPoints()->Begin() ;it != this->GetInput()->GetPoints()->End();it++)
	{

		PointType pt = it.Value();				//          std::cout <<  it.Value()<<std::endl; 
		if(this->GetOutput()->TransformPhysicalPointToIndex(it.Value(),index))
		{
			this->GetOutput()->SetPixel(index, this->GetOutput()->GetPixel(index)+1);
		}
		else
		{
			outsidePixels ++;
		}
		for(int i=0;i<3;i++)
		{
			if(pt[i]<min[i])
			{
				min[i]=pt[i];

			}
			if(pt[i]>max[i])
			{
				max[i]=pt[i];
			}

		}	
	} 
	std::cout << " ouside pixels " << outsidePixels << "min " << min << " max " << max << std::endl;

	if(outsidePixels>1000)
	{
		outsidePixels=0;	
		/*				PointType point;
						std::cout <<  " mi n" << min << std::endl;
						std::cin.get();	
						*/			 	for(int i=0;i<3;i++)
		min[i]=min[i];	

		std::cout << " min " << min << std::endl;
		std::cout << this->GetOutput()->GetOffsetTable() << min << std::endl;
		this->GetOutput()->SetOrigin(min);
		std::cin.get();
		//
		//
		/*				typename  OutputImageType::SizeType size = this->GetOutput()->GetLargestPossibleRegion().GetSize();
						IndexType start;

						for(int i=0;i<3;i++)
						start[i]=min[i];	
						typedef typename OutputImageType::RegionType RegionType;
						RegionType region(start,size);
						this->GetOutput()->SetRequestedRegion( region );*/
		std::cout << this->GetOutput()->GetRequestedRegion() << std::endl;
		this->GetOutput()->Allocate();
		this->GetOutput()->FillBuffer(0);
		for(typename MeshType::PointsContainer::Iterator it=this->GetInput()->GetPoints()->Begin() ;it != this->GetInput()->GetPoints()->End();it++)
		{
			//				std::cin.get();	
			//          std::cout <<  it.Value()<<std::endl; 
			if(this->GetOutput()->TransformPhysicalPointToIndex(it.Value(),index))
			{
				/*		for(int i = 0 ; i<3;i++)
						{
						index[i]= index[i]-min[i];
						}*/
				this->GetOutput()->SetPixel(index, this->GetOutput()->GetPixel(index)+1);
			}
			else
			{
				std::cout << index << " " <<it.Value() << std::endl;
				outsidePixels++;
			}	
		} 
		std::cout << " ouside pixels " << outsidePixels << std::endl;

	}
}

template<class TInputMesh, class TOutputImage>
float
	MeshToImageFilter<TInputMesh,TOutputImage>
::BinaryImageOfLabels(int label, int flip)
{

	float averageX=0;
	this->m_usingLabels =true;
	int count =0;
	typename MeshType::CellDataContainer::ConstIterator cellData;
	//if(label != -1 )
	cellData= this->GetInput()->GetCellData()->Begin();
	typename MeshType::PointType pt;
	IndexType index;
	for(typename MeshType::CellsContainer::ConstIterator cells = this->GetInput()->GetCells()->Begin();cells!= this->GetInput()->GetCells()->End(); cells++,++cellData)
	{
		if(label ==-1 || cellData.Value() == label )
		{
			for(typename MeshType::CellTraits::PointIdIterator  pointIdIt  =cells.Value()->PointIdsBegin();pointIdIt != cells.Value()->PointIdsEnd(); pointIdIt++)
			{

				this->GetInput()->GetPoint (*pointIdIt, &pt);
				if(this->GetOutput()->TransformPhysicalPointToIndex(pt,index))
				{
					//this->GetOutput()->SetPixel(index, 1);
					if ( flip != 0)
					{
						index[0] = flip*2 - index[0];
						this->GetOutput()->SetPixel(index, 1);
					}	
					else
					{
						this->GetOutput()->SetPixel(index, 1);

					}	
					averageX+= index[0];
					count++;
				}else
				{
					std::cout << " outside pixel " << std::endl;
				}
			} 
		}
	}
	return averageX/count;
	//std::cout << " count " <<count <<std::endl;
}



#endif
