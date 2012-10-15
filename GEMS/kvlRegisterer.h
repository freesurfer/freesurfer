#ifndef __kvlRegisterer_h
#define __kvlRegisterer_h

#include "itkImage.h"
#include "itkAffine3DTransform.h"
#include "itkShrinkImageFilter.h"
#include "itkRecursiveMultiResolutionPyramidImageFilter.h"
#include "kvlParameterOrderPowellOptimizer.h"
#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkFixedArray.h"


namespace kvl
{



/**
 *
 */
class Registerer: public itk::Object
{
public :

  /** Standard class typedefs */
  typedef Registerer  Self;
  typedef itk::Object  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( Registerer, itk::Object );

  /** Some typedefs */
  typedef itk::Image< short, 3 >  ImageType;
  typedef itk::Affine3DTransform< double >   TransformType;
  typedef TransformType::ParametersType  ParametersType;
  typedef itk::Image< float, 3 >  InternalImageType;
  typedef itk::ShrinkImageFilter< ImageType, InternalImageType >   ImageShrinkerType;
  typedef itk::RecursiveMultiResolutionPyramidImageFilter< InternalImageType,
          InternalImageType >   ImagePyramidType;
  typedef itk::MultiResolutionImageRegistrationMethod< InternalImageType,
          InternalImageType >   RegistrationType;

  // Set/get fixed image
  itkSetConstObjectMacro( FixedImage, ImageType );
  itkGetConstObjectMacro( FixedImage, ImageType );

  // Set/get moving image
  itkSetConstObjectMacro( MovingImage, ImageType );
  itkGetConstObjectMacro( MovingImage, ImageType );

  //
  void  StartRegistration();

  //
  void  ApplyParameters( ImageType* image ) const
  {
    std::vector< ImageType::Pointer >  images;
    images.push_back( image );
    this->ApplyParameters( images );
  }

  void  ApplyParameters( std::vector< ImageType::Pointer > images ) const;

  // Set/Get
  void  SetParameters( const ParametersType& parameters )
  {
    m_Transform->SetParameters( parameters );
  }

  const ParametersType&  GetParameters() const
  {
    return m_Transform->GetParameters();
  }


  // Standard Set/Get access to data members
  itkGetConstObjectMacro( FixedShrinker, ImageShrinkerType );
  itkGetConstObjectMacro( MovingShrinker, ImageShrinkerType );

  itkGetConstObjectMacro( FixedImagePyramid, ImagePyramidType );
  itkGetConstObjectMacro( MovingImagePyramid, ImagePyramidType );

  itkGetConstObjectMacro( Optimizer, ParameterOrderPowellOptimizer );

  itkGetConstObjectMacro( Registration, RegistrationType );

  itkSetMacro( UseDefaultSchedule, bool );
  itkGetConstReferenceMacro( UseDefaultSchedule, bool );

  itkSetMacro( Reseed, bool );
  itkGetConstReferenceMacro( Reseed, bool );

  itkSetMacro( NumberOfBins, int );
  itkGetConstReferenceMacro( NumberOfBins, int );

  itkSetMacro( NumberOfSamples, int );
  itkGetConstReferenceMacro( NumberOfSamples, int );

  itkSetMacro( DegreesOfFreedom, int );
  itkGetConstReferenceMacro( DegreesOfFreedom, int );

  //itkSetMacro( AcceptNewDirections, bool );
  //itkGetConstReferenceMacro( AcceptNewDirections, bool );

  itkSetMacro( MaximumNumberOfIterations, int );
  itkGetConstReferenceMacro( MaximumNumberOfIterations, int );

  itkSetMacro( InitialBracketStepSize, float );
  itkGetConstReferenceMacro( InitialBracketStepSize, float );

  //itkSetMacro( InitialBracketStepSizeShrinkFactor, float );
  //itkGetConstReferenceMacro( InitialBracketStepSizeShrinkFactor, float );

  //itkSetMacro( MaximumBracketStep, float );
  //itkGetConstReferenceMacro( MaximumBracketStep, float );

  itkSetMacro( AbsolutePrecisionBrent, float );
  itkGetConstReferenceMacro( AbsolutePrecisionBrent, float );

  itkSetMacro( MaximumNumberOfIterationsBrent, int );
  itkGetConstReferenceMacro( MaximumNumberOfIterationsBrent, int );

  //itkSetMacro( AbsolutePrecision, float );
  //itkGetConstReferenceMacro( AbsolutePrecision, float );


protected:
  Registerer();
  virtual ~Registerer();


private:
  Registerer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  // Basis data members
  ImageType::ConstPointer  m_FixedImage;
  ImageType::ConstPointer  m_MovingImage;

  TransformType::Pointer   m_Transform;


  // Data members whos address can be returned so user can observe them
  ImageShrinkerType::Pointer  m_FixedShrinker;
  ImageShrinkerType::Pointer  m_MovingShrinker;

  ImagePyramidType::Pointer  m_FixedImagePyramid;
  ImagePyramidType::Pointer  m_MovingImagePyramid;

  ParameterOrderPowellOptimizer::Pointer  m_Optimizer;

  RegistrationType::Pointer m_Registration;


  // Some convenience methods
  typedef itk::FixedArray< unsigned int, 3 >  ShrinkFactorsType;
  ShrinkFactorsType  GetShrinkFactors( const ImageType* image ) const;

  typedef ImagePyramidType::ScheduleType  ScheduleType;
  ScheduleType  GetSchedule( const ImageType* image ) const;

  typedef ParameterOrderPowellOptimizer::ParameterOrderType  ParameterOrderType;
  ParameterOrderType  GetParameterOrder() const;


  // Data members regarding loads of interal optimization params
  bool  m_UseDefaultSchedule;
  bool  m_Reseed;
  int  m_NumberOfBins;
  int  m_NumberOfSamples;
  int  m_DegreesOfFreedom;
  //bool  m_AcceptNewDirections;
  int  m_MaximumNumberOfIterations;
  float  m_InitialBracketStepSize;
  //float  m_InitialBracketStepSizeShrinkFactor;
  //float  m_MaximumBracketStep;
  float  m_AbsolutePrecisionBrent;
  int  m_MaximumNumberOfIterationsBrent;
  //float  m_AbsolutePrecision;


};


} // end namespace kvl

#endif
