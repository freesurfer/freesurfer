
// STL
#include <cmath>

// ITK
#include <itkPointSetToImageFilter.h>
#include <itkNeighborhoodConnectedImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

// VTK
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>

#include "Preprocessor.h"

template<class T>
T sqr(T x)
{ 
  return x*x; 
}

Preprocessor::Preprocessor()
{
  inputData = vtkPolyData::New();
  paddingSize.Fill(5);
  imageResolution.Fill(1.0);
  mask = NULL;
  data = NULL;

  zeroLevelValue = 100;
  insideValue = 200;
  outsideValue = 0;
}

void Preprocessor::SetPadding(const SizeType& padding)
{
  paddingSize = padding;
}

void Preprocessor::SetResolution(const ResolutionType& resolution)
{
  imageResolution = resolution;
}

void Preprocessor::SetInput(vtkPolyData* polyData)
{
  inputData->DeepCopy(polyData);
}

Preprocessor::PointSetPointer
Preprocessor::convertInputToPointSet()
{
  typedef PointSetType::PointType PointType;
  PointSetPointer pointSet = PointSetType::New();

  PointType pt;
  unsigned int ptSetCounter(0);

  // recover data by cells
  // cell 0 = 0 line
  // cell 1 = 1 line
  // cell 2 = 0 -> 1
  // cell 3 = 0 -> 1
  vtkCellArray* lines = inputData->GetLines();
  lines->InitTraversal();
  vtkIdType* vtkIds = NULL;
  vtkIdType  numPoints;
  unsigned int counter = 0;

  vtkFloatArray* floatArray = vtkFloatArray::New();
  floatArray->SetNumberOfComponents(1);
  floatArray->SetNumberOfTuples( inputData->GetPoints()->GetNumberOfPoints() );
  floatArray->SetName("distance");

  // get the constant potential lines
  float fBuf;
  while ( counter < 2 && lines->GetNextCell(numPoints, vtkIds) )
  {
    vtkIdType* idPtr = vtkIds;
    for (vtkIdType ui=0; ui < numPoints; ++ui, ++idPtr )
    {
      for(unsigned int uiDim = 0; uiDim < Dimension; ++uiDim)
        pt[uiDim] = inputData->GetPoints()->GetPoint(*idPtr)[uiDim];
      pointSet->SetPoint(ptSetCounter++, pt);
      fBuf = (float)counter;
      floatArray->SetTupleValue( *idPtr, &fBuf);
    } // next ui
    ++counter;
  }
  if ( counter < 2 )
  {
    std::cerr << " incorect number of cells\n";
    exit(1);
  }

  // get the transverse lines
  counter = 0;
  fBuf = .0f;
  while ( counter < 2 && lines->GetNextCell(numPoints, vtkIds) )
  {
    typedef std::map<unsigned int,float> MapType;
    MapType distMap;
    float fLength(.0f);
    vtkIdType* idPtr = vtkIds;

    // the following is useful for length measurements
    PointType prevPt;
    for(unsigned int uiDim(0); uiDim < Dimension; ++uiDim)
      prevPt[uiDim] = inputData->GetPoints()->GetPoint(*idPtr)[uiDim];

    for(vtkIdType ui=0; ui < numPoints; ++ui, ++idPtr)
    {
      for(unsigned int uiDim = 0; uiDim < Dimension; ++uiDim)
        pt[uiDim] = inputData->GetPoints()->GetPoint(*idPtr)[uiDim];
      pointSet->SetPoint(ptSetCounter,pt);

      // update length
      fBuf = .0f;
      for(unsigned int uiDim=0; uiDim < Dimension; ++uiDim)
        fBuf += sqr( pt[uiDim] - prevPt[uiDim] );
      fLength += std::sqrt( fBuf );
      distMap[ptSetCounter++] = fLength;
      prevPt = pt; // compute the distance - at each moment, remember the last point to speed up things
    } // next ui, idPtr

    for(MapType::const_iterator cit = distMap.begin();
        cit != distMap.end(); ++cit )
    {
      fBuf = cit->second / fLength;
      floatArray->SetTupleValue( cit->first, &fBuf );
    }

    counter++;
  }
  if ( counter < 2 )
  {
    std::cerr << " incorect number of transverse cells\n";
    exit(1);
  }

  inputData->GetPointData()->SetScalars( floatArray );

  return pointSet;
}

Preprocessor::LabelImagePointer
Preprocessor::convertPointSetToBinaryImage(PointSetPointer pointSet)
{
  typedef itk::PointSetToImageFilter<PointSetType, LabelImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();

  // compute the bounding box of the pointSet
  PointSetType::BoundingBoxType::BoundsArrayType bbox
    = pointSet->GetBoundingBox()->GetBounds();

  // use it to infer the resolution
  // also use the padding size
  FilterType::SizeType imageSize;
  for(unsigned int uiDim=0; uiDim < Dimension; ++uiDim)
    imageSize[uiDim] = (int)std::floor( (bbox[2*uiDim+1] - bbox[2*uiDim]) 
           / imageResolution[uiDim] ) + 2*paddingSize[uiDim]; 
  filter->SetSize(imageSize);
  
  // set the origin accordingly
  FilterType::PointType origin;
  for(unsigned int uiDim = 0; uiDim < Dimension; ++uiDim )
    origin[uiDim] = bbox[2*uiDim] - imageResolution[uiDim] * paddingSize[uiDim];
  filter->SetOrigin(origin);

  // spacing
  FilterType::SpacingType spacing;
  for(unsigned int uiDim = 0; uiDim < Dimension; ++uiDim )
    spacing[uiDim] = imageResolution[uiDim];
  filter->SetSpacing(spacing);

  // set non-zero value for pixel trace
  filter->SetInsideValue(zeroLevelValue);
  filter->SetInput(pointSet);
  filter->Update();

  return filter->GetOutput();
}

void Preprocessor::generateImageMask(LabelImagePointer image)
{
  typedef 
    itk::NeighborhoodConnectedImageFilter<LabelImageType, LabelImageType> 
    FilterType;
  FilterType::Pointer filter = FilterType::New();
  
  filter->SetInput( image );
  // set a seed
  FilterType::IndexType seed;
  seed.Fill(0);
  filter->SetSeed(seed);

  // control size of the neighborhood
  FilterType::InputImageSizeType radius;
  radius.Fill(0);
  filter->SetRadius(radius);
  
  // set threshold value
  filter->SetUpper( (int)zeroLevelValue - 1 );

  // set outside value
  filter->SetReplaceValue( 1 );

  filter->Update();

  // now form the final mask image
  typedef itk::ImageRegionIterator<LabelImageType> IteratorType;
  typedef itk::ImageRegionConstIterator<LabelImageType> ConstIteratorType;

  IteratorType iter(image, image->GetRequestedRegion() );
  ConstIteratorType citerOutside( filter->GetOutput(), 
          filter->GetOutput()->GetRequestedRegion() );
  
  for( iter.GoToBegin(), citerOutside.GoToBegin();
       !iter.IsAtEnd(); ++iter, ++citerOutside)
  {
    if ( citerOutside.Get() ) iter.Set( outsideValue );
    else if ( !iter.Get() ) iter.Set(insideValue);
  } // next iter, citerOutside

}


/*
Create a data image
On pixels corresponding to the zero level,
pass over the values in the closest point from the pointSet
*/
void Preprocessor::generateInitialDataImage()
{
  data = DataImageType::New();
  data->SetRegions( mask->GetRequestedRegion() );
  data->Allocate();
  data->CopyInformation( mask );

  data->FillBuffer(.0f);

  typedef itk::ImageRegionConstIterator<LabelImageType> ConstIteratorType;
  ConstIteratorType citer(mask, mask->GetRequestedRegion() );
  typedef itk::ImageRegionIteratorWithIndex<DataImageType> DataIteratorType;
  DataIteratorType iterData(data, data->GetRequestedRegion() );

  double pt[3];
  vtkIdType foundId;
  PointSetType::PointType point;
  for( iterData.GoToBegin(), citer.GoToBegin();
       !iterData.IsAtEnd(); ++iterData, ++citer)
  {
    if ( citer.Get() == zeroLevelValue )
    {
      data->TransformIndexToPhysicalPoint( iterData.GetIndex(), point );
      for(unsigned int uiDim = 0; uiDim < Dimension; ++uiDim )
        pt[uiDim] = point[uiDim];
      foundId = inputData->FindPoint(pt);
      // set the value of the pixel
      iterData.Set( inputData->GetPointData()->GetScalars()->GetComponent(foundId,0) );
    }
  } // next iterData, citer
}

void Preprocessor::Update()
{
  if ( !inputData )
  {
    std::cerr << " Preprocessor::Update - no input data\n";
    exit(1);
  }

  PointSetPointer pointSet = convertInputToPointSet(); // also adds a float array to the input
  mask = convertPointSetToBinaryImage(pointSet);
  // next is the in-place completion of the binContour mask
  generateImageMask(mask); // here the mask gets allocated and stuff
  generateInitialDataImage(); // uses the inputData and sets the data member
}

