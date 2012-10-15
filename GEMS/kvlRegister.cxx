#include "kvlRegisterer.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMGHImageIOFactory.h"
#include "kvlProgressReporter.h"
#include "itkCenteredTransformInitializer.h"


namespace kvl
{

class OptimizerWatcher
{
public:

  typedef itk::SimpleMemberCommand< OptimizerWatcher > CommandType;

  // Constructor
  OptimizerWatcher( const ParameterOrderPowellOptimizer* optimizer,
                    const Registerer::RegistrationType* registration )
  {
    m_Optimizer = optimizer;
    m_Registration = registration;

    m_OptimizerIterationCommand = CommandType::New();
    m_OptimizerIterationCommand->SetCallbackFunction( this,
        &OptimizerWatcher::HandleOptimizerIteration );
    m_Optimizer->AddObserver( itk::IterationEvent(), m_OptimizerIterationCommand );

    m_OptimizerEndCommand = CommandType::New();
    m_OptimizerEndCommand->SetCallbackFunction( this,
        &OptimizerWatcher::HandleOptimizerEnd );
    m_Optimizer->AddObserver( itk::EndEvent(), m_OptimizerEndCommand );


    m_RegistrationIterationCommand = CommandType::New();
    m_RegistrationIterationCommand->SetCallbackFunction( this,
        &OptimizerWatcher::HandleRegistrationIteration );
    m_Registration->AddObserver( itk::IterationEvent(), m_RegistrationIterationCommand );


  }

  //
  virtual ~OptimizerWatcher() {};

protected:

  //
  virtual void HandleOptimizerIteration()
  {
    std::cout << "       Iteration " << m_Optimizer->GetCurrentIteration()
              << ": value = " << m_Optimizer->GetCurrentCost() << std::endl;
    const ParameterOrderPowellOptimizer::ParametersType&  currentPosition = m_Optimizer->GetCurrentPosition();
    std::cout << "            Current position: " << std::endl;
    std::cout << "                    Rotation: ";
    for ( int i = 0; i < 3; i++ )
    {
      std::cout << currentPosition[ i ] * 180.0 / 3.14 << "  ";
    }
    std::cout << std::endl;
    std::cout << "                 Translation: ";
    for ( int i = 3; i < 6; i++ )
    {
      std::cout << currentPosition[ i ] << "  ";
    }
    std::cout << std::endl;
    std::cout << "                     Scaling: ";
    for ( int i = 6; i < 9; i++ )
    {
      std::cout << currentPosition[ i ] << "  ";
    }
    std::cout << std::endl;
    std::cout << "                     Skewing: ";
    for ( int i = 9; i < 12; i++ )
    {
      std::cout << currentPosition[ i ] << "  ";
    }
    std::cout << std::endl;

  }

  //
  virtual void HandleOptimizerEnd()
  {
    std::cout << m_Optimizer->GetStopConditionDescription() << std::endl;
  }

  //
  virtual void HandleRegistrationIteration()
  {
    std::cout << "Doing registration at level " << m_Registration->GetCurrentLevel()
              << " of " << m_Registration->GetNumberOfLevels() << std::endl;
  }



private:
  ParameterOrderPowellOptimizer::ConstPointer m_Optimizer;
  Registerer::RegistrationType::ConstPointer  m_Registration;

  CommandType::Pointer  m_OptimizerIterationCommand ;
  CommandType::Pointer  m_OptimizerEndCommand ;
  CommandType::Pointer  m_RegistrationIterationCommand ;

};



};




int main( int argc, char* argv[] )
{
  // Sanity check on input
  if ( argc < 3 )
  {
    std::cerr << argv[0] << " inputImage referenceImage [ degreesOfFreedom=6 initialTranslationMode=0 "
              << "numberOfBins=20 useDefaultSchedule=1 extraInputImage1 ... ]" << std::endl;
    return -1;
  }

  // Parse input
  const std::string  fixedImageFileName = argv[ 1 ];
  const std::string  movingImageFileName = argv[ 2 ];

  int  degreesOfFreedom = 6;
  if ( argc > 3 )
  {
    std::istringstream  degreesOfFreedomStream( argv[ 3 ] );
    degreesOfFreedomStream >> degreesOfFreedom;
  }
  int  initialTranslationMode = 0;
  if ( argc > 4 )
  {
    std::istringstream  initialTranslationModeStream( argv[ 4 ] );
    initialTranslationModeStream >> initialTranslationMode;
  }
  int  numberOfBins = 20;
  if ( argc > 5 )
  {
    std::istringstream  numberOfBinsStream( argv[ 5 ] );
    numberOfBinsStream >> numberOfBins;
  }
  bool  useDefaultSchedule = true;
  if ( argc > 6 )
  {
    std::istringstream  useDefaultScheduleStream( argv[ 6 ] );
    useDefaultScheduleStream >> useDefaultSchedule;
  }
  std::vector< std::string >  extraInputImageFileNames;
  if ( argc > 7 )
  {
    for ( int i = 7; i < argc; i++ )
    {
      extraInputImageFileNames.push_back( std::string( argv[ 7 ] ) );
    }
  }


  // Print out what we have
  std::cout << "fixedImageFileName: " << fixedImageFileName << std::endl;
  std::cout << "movingImageFileName: " << movingImageFileName << std::endl;
  std::cout << "degreesOfFreedom: " << degreesOfFreedom << std::endl;
  std::cout << "numberOfBins: " << numberOfBins << std::endl;
  std::cout << "useDefaultSchedule: " << useDefaultSchedule << std::endl;
  std::cout << "initialTranslationMode: " << initialTranslationMode << std::endl;
  std::cout << "extraInputImageFileNames: " << std::endl;
  for ( std::vector< std::string >::const_iterator  it = extraInputImageFileNames.begin();
        it != extraInputImageFileNames.end(); ++it )
  {
    std::cout << "    " << *it << std::endl;
  }

  // Add support for MGH file format to ITK. An alternative way to add this by default would be
  // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
  itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );

  try
  {
    // Read images
    typedef kvl::Registerer::ImageType  ImageType;
    typedef itk::ImageFileReader< ImageType >  ReaderType;

    ReaderType::Pointer  fixedReader = ReaderType::New();
    fixedReader->SetFileName( fixedImageFileName.c_str() );
    fixedReader->Update();
    ImageType::Pointer  fixedImage = fixedReader->GetOutput();

    ReaderType::Pointer  movingReader = ReaderType::New();
    movingReader->SetFileName( movingImageFileName.c_str() );
    movingReader->Update();
    ImageType::ConstPointer  movingImage = movingReader->GetOutput();

    std::vector< ImageType::Pointer >  imagesToRegister;
    imagesToRegister.push_back( fixedImage );
    std::vector< std::string >  imagesToRegisterFileNames;
    imagesToRegisterFileNames.push_back( fixedImageFileName );
    for ( std::vector< std::string >::const_iterator  it = extraInputImageFileNames.begin();
          it != extraInputImageFileNames.end(); ++it )
    {
      ReaderType::Pointer  reader = ReaderType::New();
      reader->SetFileName( it->c_str() );
      reader->Update();
      imagesToRegister.push_back( reader->GetOutput() );

      imagesToRegisterFileNames.push_back( *it );
    }

    // Set up the registerer
    kvl::Registerer::Pointer  registerer = kvl::Registerer::New();
    registerer->SetFixedImage( fixedImage );
    registerer->SetMovingImage( movingImage );
    registerer->SetDegreesOfFreedom( degreesOfFreedom );
    registerer->SetNumberOfBins( numberOfBins );
    registerer->SetUseDefaultSchedule( useDefaultSchedule );


    // Add some observers to what it will be doing
    kvl::ProgressReporter  reporter1( registerer->GetFixedImagePyramid(), "Calculating fixed image pyramid" );
    kvl::ProgressReporter  reporter2( registerer->GetMovingImagePyramid(), "Calculating moving image pyramid" );
    kvl::ProgressReporter  reporter3( registerer->GetFixedShrinker(), "Shrinking fixed image" );
    kvl::ProgressReporter  reporter4( registerer->GetMovingShrinker(), "Shrinking moving image" );
    kvl::OptimizerWatcher  optimizerWatcher( registerer->GetOptimizer(),
        registerer->GetRegistration() );


    // Initialize the gross translation component, if so desired
#if 0
    kvl::Registerer::ParametersType  parameters = registerer->GetParameters();
    parameters[ 3 ] += 30;
    parameters[ 4 ] += 40;
    parameters[ 5 ] += 50;
    registerer->SetParameters( parameters );
#endif
    if ( initialTranslationMode )
    {
      // Set up the initializer
      typedef itk::AffineTransform< double, 3 >  AffineTransformType;
      typedef itk::CenteredTransformInitializer< AffineTransformType, ImageType, ImageType >  InitializerType;

      AffineTransformType::Pointer  affineTransform = AffineTransformType::New();
      InitializerType::Pointer  initializer = InitializerType::New();
      initializer->SetTransform( affineTransform );
      initializer->SetFixedImage( fixedImage );
      initializer->SetMovingImage( movingImage );

      if ( initialTranslationMode == 1 )
      {
        initializer->GeometryOn();
      }
      else
      {
        initializer->MomentsOn();
      }

      // Let the beast go
      initializer->InitializeTransform();

      // Convert what we've found into a usable form
      typedef kvl::Registerer::TransformType  TransformType;
      TransformType::Pointer  transform = TransformType::New();
      transform->SetMatrix( affineTransform->GetMatrix() );
      transform->SetOffset( affineTransform->GetOffset() );
      kvl::Registerer::ParametersType  parameters = transform->GetParameters();

      // Initialize the registerer accordingly
      std::cout << "Initializing the parameters as follows: " << parameters << std::endl;
      registerer->SetParameters( parameters );
    }


    // Let the beast go
    if ( degreesOfFreedom > 0 )
    {
      registerer->StartRegistration();
    }

    // Apply the found parameters and write out
    for ( unsigned int i = 0; i < imagesToRegister.size(); i++ )
    {
      const std::string  fileName = imagesToRegisterFileNames[ i ];
      ImageType::Pointer  image = imagesToRegister[ i ];

      // Apply the found parameters
      registerer->ApplyParameters( image );

      // Write out
      std::ostringstream  outputFileNameStream;
      outputFileNameStream << itksys::SystemTools::GetFilenameWithoutExtension( fileName.c_str() )
                           << "_coregistered.mgz";
      const std::string  outputFileName = outputFileNameStream.str();
      typedef itk::ImageFileWriter< ImageType >  WriterType;
      WriterType::Pointer  writer = WriterType::New();
      writer->SetInput( image );
      writer->SetFileName( outputFileName.c_str() );
      writer->Update();

      std::cout << "Wrote out " << outputFileName << std::endl;
    }

  }
  catch ( itk::ExceptionObject& e )
  {
    std::cerr << e << std::endl;
    return -1;
  }


  return 0;
};

