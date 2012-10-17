/*

Gheorghe Postelnicu, 2006

Laplace pre-processor class

Takes a vtk polydata object and returns a binary image
such that the integration goes smoothly.

Input data is assumed to be as follows:
- if 2 dimensional data is to be used, then DO FORMAT the data prior to that
 (that means, discard the useless coordinate)

- in its current form, the data is assumed to contain 4 cells, as follows:
 o the first 2 cells will be potential iso-lines;
 o the last 2 cells will be potential traversal lines;

*/
#ifndef _Preprocessor_h
#define _Preprocessor_h

#include <vtkPolyData.h>
#include <itkImage.h>
#include <itkPointSet.h>
#include <itkFixedArray.h>

class Preprocessor
{
 public:
  static const unsigned int Dimension = 2;
  typedef unsigned char LabelPixelType;
  typedef double DataPixelType;

  typedef itk::Image<LabelPixelType, Dimension> LabelImageType;
  typedef LabelImageType::Pointer LabelImagePointer;
  typedef itk::Image<DataPixelType, Dimension>  DataImageType;
  typedef DataImageType::Pointer DataImagePointer;

  typedef itk::FixedArray<unsigned int,Dimension> SizeType;
  typedef itk::FixedArray<float, Dimension> ResolutionType;

  Preprocessor();
  void SetInput(vtkPolyData* data);
  void SetPadding(const SizeType&);
  void SetResolution(const ResolutionType&);
  void Update();
  inline LabelImagePointer GetOutputLabel() { return mask; }
  inline DataImagePointer  GetOutputData() { return data; }

  inline LabelPixelType GetMaskInsideValue()
    { return insideValue; }
  inline LabelPixelType GetMaskZeroValue()
    { return zeroLevelValue; }

 private:
  typedef itk::PointSet<float,Dimension> PointSetType;
  typedef PointSetType::Pointer PointSetPointer;
  PointSetPointer convertInputToPointSet();
  
  LabelImagePointer convertPointSetToBinaryImage(PointSetPointer);
  void              generateImageMask(LabelImagePointer);
  void              generateInitialDataImage();

  vtkPolyData*    inputData;
  SizeType        paddingSize;
  ResolutionType  imageResolution;

  LabelImagePointer mask;
  DataImagePointer  data;
  
  LabelPixelType zeroLevelValue;
  LabelPixelType insideValue;
  LabelPixelType outsideValue;
};

#endif
