
#include <iostream>

#define export // obsolete feature "export template" used in these header files
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>
#include <itkDerivativeImageFilter.h>
#undef export

#include "Tracer.h"

bool Tracer::DoNotExitOnError = false;

Tracer::Tracer()
{
  mask = NULL;
  data = NULL;

  imageInterpolator = NULL;
  for(unsigned int ui=0; ui<Dimension; ++ui)
    derivativeInterpolator[ui] = NULL;
}

Tracer::~Tracer()
{
}

/*
set input image
- also compute the interpolators
*/
void
Tracer::SetInputData(ImagePointer image)
{
  data = ImageSOType::New();
  data->SetImage(image);

  // set the stepSize
  ImageType::SpacingType spacing
    = image->GetSpacing();
  stepSize = std::min( spacing[0], spacing[1] ) * 2;

  // set interpolators
  typedef itk::DerivativeImageFilter<ImageType,ImageType> FilterType;
  FilterType::Pointer filter[Dimension];

  for(unsigned int ui=0; 
      ui < Dimension; ++ui)
    {
      filter[ui] = FilterType::New();
      filter[ui]->SetDirection(ui);
      filter[ui]->SetInput( data->GetImage() );
      filter[ui]->Update();

      derivativeInterpolator[ui] = InterpolatorType::New();
      derivativeInterpolator[ui]->SetInputImage( filter[ui]->GetOutput() );
    } // next ui

  imageInterpolator = InterpolatorType::New();
  imageInterpolator->SetInputImage( data->GetImage() );
}

void
Tracer::SetInputMask(MaskImagePointer inputMask)
{
  mask = MaskImageSOType::New();
  mask->SetImage( inputMask );
}

/** Creates a vtkPointLocator object from the polydata points */
void
Tracer::SetInputContours(vtkPolyData* data)
{
  std::cout << "Set input contours\n";
  pts = data->GetPoints();
  locator = vtkPointLocator::New();

  locator->SetDataSet(data);

  locator->BuildLocator();
  locator->AutomaticOn();
  std::cout << "\t DONE\n";
}


Tracer::LineType
Tracer::ComputeMidline() const
{
  return this->ComputeIsoline(.5);
}

Tracer::LineType
Tracer::ComputeIsoline(double dval) const
{
  // start by recovering a seed
  PointType seed = this->GetClosestPoint(dval);
  LineType line;
  line.push_back(seed);
  
  bool activeBegin(true), activeEnd(true);

  PointType ptBuf;
  while( activeBegin || activeEnd )
  {
    if (activeBegin)
    {
      activeBegin = this->StepTangent(line.front(), ptBuf, dval);
      if (activeBegin) line.push_front(ptBuf);
    }
    if (activeEnd)
    {
      activeEnd = this->StepTangent(line.back(), ptBuf, dval, true);
      if (activeEnd) line.push_back(ptBuf);
    }
  } 

  // compute and print stats
  double dmin(1000), dmax(-1), dsum(0), dbuf;
  for( LineType::const_iterator cit = line.begin();
       cit != line.end(); ++cit )
  {
    this->ValueAt(*cit, dbuf);
    dmin = std::min(dmin, dbuf);
    dmax = std::max(dmax, dbuf);
    dsum += dbuf;
  }
  dsum /= line.size();
  std::cout << " target = " << dval << std::endl;
  std::cout << " min    = " << dmin << std::endl;
  std::cout << " max    = " << dmax << std::endl;
  std::cout << " avg    = " << dsum << std::endl;

  return line;

}

/*

the orientation is important for the profiles,
since they are the actual goal of all this

the positive order should be in the increasing order of the gradient

*/
Tracer::LineType
Tracer::ComputeProfile(double x, double y) const
{
  PointType seed;
  seed[0] = x;
  seed[1] = y;

  LineType line;
  line.push_back(seed);

  bool activeBegin(true), activeEnd(true);
  PointType ptBuf;
  while( (activeBegin || activeEnd) && line.size()<500 )
  {
    if (activeBegin)
    {
      activeBegin = this->StepGradient( line.front(), ptBuf, true);
      line.push_front( ptBuf );
    }
    if (activeEnd)
    {
      activeEnd = this->StepGradient( line.back(), ptBuf);
      line.push_back( ptBuf );
    }
  }

  return line;
}


void
Tracer::ValueAt(const PointType& pt, double& dval) const
{
  InterpolatorType::PointType cpt;
  std::copy( pt.Begin(), pt.End(), cpt.Begin() );

  dval = imageInterpolator->Evaluate( cpt );
}

void
Tracer::DerivativeAt(const PointType& pt, VectorType& derivative) const
{
  InterpolatorType::PointType cpt;
  std::copy( pt.Begin(), pt.End(), cpt.Begin() );

  for(unsigned int ui=0;
      ui < Dimension; ++ui)
    derivative[ui] = derivativeInterpolator[ui]->Evaluate( cpt );
}

Tracer::VectorType
Tracer::GetTangent(const PointType& pt) const
{
  VectorType gradient;
  this->DerivativeAt(pt, gradient);

  // compute the perpendicular
  VectorType tangent;
  tangent[0] = - gradient[1];
  tangent[1] =   gradient[0];

  tangent.Normalize();
  tangent *= stepSize;
  
  return tangent;
}

Tracer::FrameType
Tracer::GetLocalFrame(const PointType& pt) const
{
  VectorType gradient;
  this->DerivativeAt(pt, gradient);
  gradient.Normalize();
  gradient *= stepSize;

  // compute the tangent
  VectorType tangent;
  tangent[0] = - gradient[1];
  tangent[1] =   gradient[0];

  return std::make_pair(tangent, gradient);
}

bool
Tracer::StepGradient(const PointType& pt,
         PointType& returnedPoint,
         bool negative) const
{
  // use the locator to see if close to boundary
  double buf[3];
  buf[0] = pt[0];
  buf[1] = pt[1];
  buf[2] = .0;
  double dist2;

  vtkIdType id = locator->FindClosestPointWithinRadius(1.41 *stepSize, buf, dist2);
  if (id >= 0)
  {
    for(unsigned int ui=0; ui<Dimension; ++ui)
      returnedPoint[ui] = pts->GetPoint(id)[ui];
    return false;
  }

  FrameType frame = this->GetLocalFrame(pt);
  
  //VectorType& tangent = frame.first;
  VectorType& gradient= frame.second;
  if (negative) gradient *= -1;

  for(unsigned int ui=0; ui<Dimension; ++ui)
    returnedPoint[ui] = pt[ui] + gradient[ui];

#if ITK_VERSION_MAJOR >= 5  
  bool inside = mask->IsInsideInWorldSpace( returnedPoint );
#else  
  bool inside = mask->IsInside( returnedPoint );
#endif  
  // this will connect all lines to the boundary but
  // may cause intersection
  if (! inside)
  {
     double buf[3];
     //buf[0] = pt[0];
     //buf[1] = pt[1];
     buf[0] = returnedPoint[0];
     buf[1] = returnedPoint[1];
     buf[2] = .0;
     //double dist2;
 
     vtkIdType id = locator->FindClosestPoint(buf);
     if (id >= 0)
     {
       for(unsigned int ui=0; ui<Dimension; ++ui)
         returnedPoint[ui] = pts->GetPoint(id)[ui];
       // return false;
     }
  
  }

  return inside;
}


bool
Tracer::StepTangent(const PointType& pt,
        PointType& returnedPoint,
        double dval,
        bool negative) const
{
  FrameType frame = this->GetLocalFrame(pt);

  VectorType& tangent = frame.first;
  VectorType& gradient = frame.second;
  if (negative) tangent *= -1;

  for(unsigned int ui=0; ui<Dimension; ++ui)
    returnedPoint[ui] = pt[ui] + tangent[ui];

#if ITK_VERSION_MAJOR >= 5  
  if (!mask->IsInsideInWorldSpace(returnedPoint)) return false;
#else
  if (!mask->IsInside(returnedPoint)) return false;
#endif  

  // corect the trajectory - use the gradient and a first order aproximation
  double dbuf;
  ImageSOType::PointType ptBuf = returnedPoint;
  this->ValueAt(ptBuf, dbuf);
  double dif = dbuf - dval;
  for(unsigned int ui=0; ui<Dimension; ++ui)
    returnedPoint[ui] -= gradient[ui] * dif;

  this->ValueAt(ptBuf,dbuf);
  if ( std::abs(dbuf - dval) < dif ) returnedPoint = ptBuf;

#if ITK_VERSION_MAJOR >= 5  
  return mask->IsInsideInWorldSpace( returnedPoint );
#else  
  return mask->IsInside( returnedPoint );
#endif  
}

Tracer::PointType
Tracer::GetClosestPoint(double val) const
{
  typedef itk::ImageRegionConstIteratorWithIndex<ImageType> ImageConstIteratorWithIndex;
  typedef itk::ImageRegionConstIterator<MaskImageType> ImageConstIterator;
  ImageConstIteratorWithIndex cit( data->GetImage(), data->GetImage()->GetRequestedRegion() );
  ImageConstIterator citMask( mask->GetImage(), mask->GetImage()->GetRequestedRegion() );

  cit.GoToBegin();
  ImageType::IndexType argMin(cit.GetIndex());
  double dMin = 10.0;

  for(cit.GoToBegin(), citMask.GoToBegin(); !cit.IsAtEnd(); ++cit, ++citMask)
    {
      if ( citMask.Get() )
  if ( std::abs( val - cit.Get() ) < dMin )
    {
      dMin = std::abs( val - cit.Get() );
      argMin = cit.GetIndex();
    }
    } // next cit, citMask

  if ( dMin > 9.0 )
    {
      std::cerr << " error finding close point while searching for value " << val << "\n";
      if (!DoNotExitOnError) exit(1);
    }

  PointType retVal;
  data->GetImage()->TransformIndexToPhysicalPoint(argMin, retVal);
  return retVal;
}
