#include "kvlRegisterer.h"

#include "itkAffineTransform.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "vnl/vnl_sample.h"


namespace kvl
{

//
//
//
Registerer
::Registerer()
{
  m_FixedImage = 0;
  m_MovingImage = 0;

  m_Transform = TransformType::New();

  m_UseDefaultSchedule = true;

  m_FixedShrinker = ImageShrinkerType::New();
  m_MovingShrinker = ImageShrinkerType::New();

  m_FixedImagePyramid = ImagePyramidType::New();
  m_MovingImagePyramid = ImagePyramidType::New();

  m_Optimizer  = ParameterOrderPowellOptimizer::New();

  m_Registration = RegistrationType::New();

  // Default values for loads of internal params
  m_Reseed = false;
  m_NumberOfBins = 20;
  m_NumberOfSamples = 20000;
  m_DegreesOfFreedom = 6;
  //m_AcceptNewDirections = false;
  m_MaximumNumberOfIterations = 10;
  m_InitialBracketStepSize = 5.0f;
  //m_InitialBracketStepSizeShrinkFactor = 2.0f;
  //m_MaximumBracketStep = 30.0f;
  m_AbsolutePrecisionBrent = 0.1;
  m_MaximumNumberOfIterationsBrent = 50;
  //m_AbsolutePrecision = 0.1;

}




//
//
//
Registerer
::~Registerer()
{

}



//
//
//
void
Registerer
::StartRegistration()
{

  // Sanity check on input
  if ( !m_FixedImage || !m_MovingImage )
  {
    itkExceptionMacro( << "Can't run without having both a fixed image and a moving image" );
  }


  // Shrink and simultaneously cast input images to float
  ShrinkFactorsType  fixedShrinkFactors = this->GetShrinkFactors( m_FixedImage );
  std::cout << "Using fixedShrinkFactors: " << fixedShrinkFactors << std::endl;
  m_FixedShrinker->SetInput( m_FixedImage );
  for ( int i = 0; i < 3; i++ )
  {
    m_FixedShrinker->SetShrinkFactor( 0, fixedShrinkFactors[ i ] );
  }
  m_FixedShrinker->Update();

  ShrinkFactorsType  movingShrinkFactors = this->GetShrinkFactors( m_MovingImage );
  std::cout << "Using movingShrinkFactors: " << movingShrinkFactors << std::endl;
  m_MovingShrinker->SetInput( m_MovingImage );
  for ( int i = 0; i < 3; i++ )
  {
    m_MovingShrinker->SetShrinkFactor( 0, movingShrinkFactors[ i ] );
  }
  m_MovingShrinker->Update();


  // Set up image pyramids
  const ScheduleType  fixedSchedule = this->GetSchedule( m_FixedImage );
  const int numberOfLevels = fixedSchedule.rows();
  m_FixedImagePyramid->SetNumberOfLevels( numberOfLevels );
  m_FixedImagePyramid->SetSchedule( fixedSchedule );
  std::cout << "Using fixedSchedule: " << fixedSchedule << std::endl;

  const ScheduleType  movingSchedule = this->GetSchedule( m_MovingImage );
  m_MovingImagePyramid->SetNumberOfLevels( numberOfLevels );
  m_MovingImagePyramid->SetSchedule( movingSchedule );
  std::cout << "Using movingSchedule: " << movingSchedule << std::endl;


  // Set up metric
  typedef itk::MattesMutualInformationImageToImageMetric< InternalImageType,
          InternalImageType >   MetricType;
  if ( m_Reseed )
  {
    vnl_sample_reseed( 12345 ); // Make sure we get the same samples every time
  }
  MetricType::Pointer  metric = MetricType::New();
  metric->SetNumberOfHistogramBins( m_NumberOfBins );
  metric->SetNumberOfSpatialSamples( m_NumberOfSamples );


  // Set up m_Optimizer
  ParameterOrderPowellOptimizer::ScalesType scales( 12 );
  scales[0] = 1.0 / 0.0175;
  scales[1] = 1.0 / 0.0175;
  scales[2] = 1.0 / 0.0175;
  scales[3] = 1.0 / 1.0;
  scales[4] = 1.0 / 1.0;
  scales[5] = 1.0 / 1.0;
  scales[6] = 1.0 / 0.01;
  scales[7] = 1.0 / 0.01;
  scales[8] = 1.0 / 0.01;
  scales[9] = 1.0 / 0.01;
  scales[10] = 1.0 / 0.01;
  scales[11] = 1.0 / 0.01;
  std::cout << "Using scales: " << scales << std::endl;

  ParameterOrderType  parameterOrder = this->GetParameterOrder();
  std::cout << "Using parameterOrder: " << parameterOrder << std::endl;

  m_Optimizer->SetScales( scales );
  m_Optimizer->SetParameterOrder( parameterOrder );
  //m_Optimizer->SetAcceptNewDirections( m_AcceptNewDirections );
  //m_Optimizer->SetMaximumNumberOfIterations( m_MaximumNumberOfIterations );
  m_Optimizer->SetMaximumIteration( m_MaximumNumberOfIterations );
  //m_Optimizer->SetInitialBracketStepSize( m_InitialBracketStepSize );
  m_Optimizer->SetStepLength( m_InitialBracketStepSize );
  //m_Optimizer->SetInitialBracketStepSizeShrinkFactor( m_InitialBracketStepSizeShrinkFactor );
  //m_Optimizer->SetMaximumBracketStep( m_MaximumBracketStep );
  //m_Optimizer->SetFractionalPrecisionBrent( 0.0 ); // Make sure this is never satisfied
  //m_Optimizer->SetAbsolutePrecisionBrent( m_AbsolutePrecisionBrent );
  m_Optimizer->SetStepTolerance( m_AbsolutePrecisionBrent );
  //m_Optimizer->SetMaximumNumberOfIterationsBrent( m_MaximumNumberOfIterationsBrent );
  m_Optimizer->SetMaximumLineIteration( m_MaximumNumberOfIterationsBrent );
  //m_Optimizer->SetFractionalPrecision( 0.0 ); // Make sure this is never satisfied
  m_Optimizer->SetValueTolerance( 0.0 ); // Make sure this is never satisfied
  //m_Optimizer->SetAbsolutePrecision( m_AbsolutePrecision );
  //m_Optimizer->MinimizeOn();
  m_Optimizer->MaximizeOff();



  // Set up m_Registration
  typedef itk::LinearInterpolateImageFunction< InternalImageType, double >  InterpolatorType;
  m_Registration->SetTransform( m_Transform );
  m_Registration->SetInterpolator( InterpolatorType::New() );
  m_Registration->SetMetric( metric );
  m_Registration->SetFixedImagePyramid( m_FixedImagePyramid );
  m_Registration->SetMovingImagePyramid( m_MovingImagePyramid );
  m_Registration->SetInitialTransformParameters( this->GetParameters() );
  m_Registration->SetNumberOfLevels( numberOfLevels );
  m_Registration->SetOptimizer( m_Optimizer );
  m_Registration->SetFixedImage( m_FixedShrinker->GetOutput() );
  m_Registration->SetMovingImage( m_MovingShrinker->GetOutput() );
  m_Registration->SetFixedImageRegion( m_FixedShrinker->GetOutput()->GetBufferedRegion() );


  // Let the beast go
  m_Registration->StartRegistration();


  // Retrieve the parameters
  this->SetParameters( m_Registration->GetLastTransformParameters() );

}



//
//
//
void
Registerer
::ApplyParameters( std::vector< ImageType::Pointer > images ) const
{

  //std::cout << "Applying transformation: " << m_Transform << std::endl;

  // Loop over all images
  for ( std::vector< ImageType::Pointer >::iterator  it = images.begin();
        it != images.end(); ++it )
  {
    // Get the image
    ImageType::Pointer  image = *it;

    // Get the original image-to-world tranform of the image
    typedef itk::AffineTransform< double, 3 >  AffineTransformType;
    AffineTransformType::MatrixType  scale;
    AffineTransformType::OffsetType  offset;
    for ( int i = 0; i < 3; i++ )
    {
      scale[ i ][ i ] = image->GetSpacing()[ i ];
      offset[ i ] = image->GetOrigin()[ i ];
    }
    AffineTransformType::Pointer  original = AffineTransformType::New();
    original->SetMatrix( image->GetDirection() * scale );
    original->SetOffset( offset );
    //std::cout << "Original image-to-world: " << original << std::endl;

    // Pre-multiply it with the calculated transformation
    AffineTransformType::Pointer updated = AffineTransformType::New();
    updated->SetMatrix( m_Transform->GetMatrix() );
    updated->SetOffset( m_Transform->GetOffset() );
    updated->Compose( original, true );
    //std::cout << "Updated image-to-world: " << updated << std::endl;

    // Now overwrite image information to get exactly that effect
    ImageType::PointType   newOrigin;
    ImageType::SpacingType  newSpacing;
    ImageType::DirectionType  newDirection;
    for ( int i = 0; i < 3; i++ )
    {
      // Offset part
      newOrigin[ i ] = updated->GetOffset()[ i ];

      // For every column, determine norm (which will be voxel spacing), and normalize direction
      double  normOfColumn = 0.0;
      for ( int j = 0; j < 3; j++ )
      {
        normOfColumn += pow( updated->GetMatrix()[ j ][ i ], 2 );
      }
      normOfColumn = sqrt( normOfColumn );
      newSpacing[ i ] = normOfColumn;
      for ( int j = 0; j < 3; j++ )
      {
        newDirection[ j ][ i ] = updated->GetMatrix()[ j ][ i ] / normOfColumn;
      }
    }
    image->SetOrigin( newOrigin );
    image->SetSpacing( newSpacing );
    image->SetDirection( newDirection );

  } // end loop over all images


}



//
//
//
Registerer::ShrinkFactorsType
Registerer
::GetShrinkFactors( const ImageType* image ) const
{

  ShrinkFactorsType  shrinkFactors;
  shrinkFactors.Fill( 1 );
  if ( m_UseDefaultSchedule )
  {
    for ( int i = 0; i < 3; i++ )
    {
      shrinkFactors[ i ] = vnl_math_rnd( 1.0f / ( image->GetSpacing()[ i ] ) );
      if ( shrinkFactors[ i ] < 1 )
      {
        shrinkFactors[ i ] = 1;
      }
    }
  }

  return shrinkFactors;
}



//
//
//
Registerer::ScheduleType
Registerer
::GetSchedule( const ImageType* image ) const
{

  // Check if doing default thing
  if ( !m_UseDefaultSchedule )
  {
    ScheduleType  schedule( 1, 3 );
    for ( int i = 0; i < 3; i++ )
    {
      schedule[ 0 ][ i ] = 1;
    }

    return schedule;
  }

  // Default case
  ScheduleType  schedule( 2, 3 );
  const ShrinkFactorsType  shrinkFactors =  this->GetShrinkFactors( image );
  for ( int i = 0; i < 3; i++ )
  {
    schedule[ 1 ][ i ] = 1;

    schedule[ 0 ][ i ] = vnl_math_rnd( 4.0 / ( shrinkFactors[ i ] * image->GetSpacing()[ i ] ) );
    if ( schedule[ 0 ][ i ] < 1 )
    {
      schedule[ 0 ][ i ] = 1;
    }
  }

  return schedule;


}


//
//
//
Registerer::ParameterOrderType
Registerer
::GetParameterOrder() const
{
  if ( m_DegreesOfFreedom > 12 )
  {
    itkExceptionMacro( "Maximum allowed degrees of freedom is 12" );
  }

  ParameterOrderType  parameterOrder( 12 );
  parameterOrder.Fill( 0 );
  for ( int i = 0; i < m_DegreesOfFreedom; i++ )
  {
    if ( i < 3 )
    {
      parameterOrder[ 3 + i ] = i+1;
    }
    else if ( i < 6 )
    {
      parameterOrder[ i-3 ] = i+1;
    }
    else
    {
      parameterOrder[ i ] = i+1;
    }

  }

  return parameterOrder;
}


} // end namespace kvl
