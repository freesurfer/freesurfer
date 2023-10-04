
#ifndef _Tracer_h
#define _Tracer_h

// STL
#include <algorithm>
#include <list>

// VTK
#include <vtkPolyData.h>
#include <vtkPointLocator.h>

// ITK
#define export // obsolete feature "export template" used in these header files
#include <itkImage.h>
#include <itkImageSpatialObject.h>
#include <itkImageMaskSpatialObject.h>
#include <itkLinearInterpolateImageFunction.h>
#undef export

class Tracer
{
 public:
  // some typedefs
  static const unsigned int Dimension = 2;
  typedef double PixelType;
  typedef itk::Image<PixelType, Dimension> ImageType;
  typedef ImageType::Pointer ImagePointer;
  typedef itk::ImageSpatialObject<Dimension,PixelType> ImageSOType;
  typedef ImageSOType::Pointer ImageSOPointer;

  typedef unsigned char MaskPixelType;
  typedef itk::Image<MaskPixelType,Dimension> MaskImageType;
  typedef MaskImageType::Pointer MaskImagePointer;
  typedef itk::ImageMaskSpatialObject<Dimension> MaskImageSOType;
  typedef MaskImageSOType::Pointer MaskImageSOPointer;

  typedef MaskImageSOType::PointType PointType;
  typedef std::list<PointType> LineType;

  Tracer();
  ~Tracer();

  void SetInputData(ImagePointer image);
  void SetInputMask(MaskImagePointer inputMask);
  void SetInputContours(vtkPolyData* contours);

  LineType ComputeMidline() const;
  LineType ComputeIsoline(double dval) const;
  LineType ComputeProfile(double x, double y) const;

  static bool DoNotExitOnError;
  
 private:
  typedef ImageSOType::CovariantVectorType VectorType;
  typedef std::pair<VectorType,VectorType> FrameType;

  // initial region contour - used to solve for profiles
  vtkPointLocator* locator;
  vtkPoints* pts;
  
  // returns tangent at a point
  VectorType GetTangent(const PointType& pt) const; 
  FrameType  GetLocalFrame(const PointType& pt) const;
  bool StepTangent(const PointType&,
		   PointType& returnedPoint,
		   double dval,
		   bool negative=false) const;
  bool StepGradient(const PointType&,
		    PointType& returnedPoint,
		    bool negative=false) const;
  float stepSize;
  MaskImageSOPointer mask;
  ImageSOPointer data;

  typedef itk::LinearInterpolateImageFunction<ImageType> InterpolatorType;
  typedef InterpolatorType::Pointer InterpolatorPointer;
  InterpolatorPointer imageInterpolator;
  InterpolatorPointer derivativeInterpolator[Dimension];
  
  void ValueAt(const PointType& pt, double& dval) const;
  void DerivativeAt(const PointType& pt, VectorType& derivative) const;

  // iterate through the image and get the closest point to a certain value
  PointType GetClosestPoint(double) const;
};

#endif
