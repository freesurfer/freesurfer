#include <iostream>
#include <string>
#include <cstdlib>
#include <unistd.h>
#include <sys/utsname.h>

// for reading in the cosine directions because itk doesn't do it properly
#include <nifti1_io.h>

#include <itkImageSeriesReader.h>
#include <itkImageFileWriter.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkNiftiImageIO.h>
#include <itkOrientedImage.h>

#include "io/itkMGHImageIOFactory.h"
#include "io/SymmetricTensorReaderStrategy.h"
#include "io/AsymmetricTensorReaderStrategy.h"
#include "io/AsymmetricTensorVectorReaderStrategy.h"

#include "datamodel/utils/itkPoistatsFilter.h"

#include "datamodel/events/CommandUpdate.h"

#include "PoistatsCLICLP.h"

bool
IsNifti( std::string fileExtension ) {

  bool isNifti = false;
  
  if( fileExtension == ".nii" || fileExtension == ".nii.gz" ) {
    isNifti = true;
  }
  
  return isNifti;
}

bool
IsMGH( std::string fileExtension ) {

  bool isMGH = false;
  
  if( fileExtension == ".mgh" || fileExtension == ".mgz" ) {
    isMGH = true;
  }
  
  return isMGH;
}

mat44
GetNiftiTransform( const char *filename ) {
  mat44 rval;
  nifti_image *img = nifti_image_read(filename, false);
  if(img == 0) {
    return rval;
  }
  
  // TODO: check the transform
  rval = img->qto_xyz; // or img->sto_xyz -- you need to check the transform type
                                     // to ensure getting the 'one true transform' for this file
  nifti_image_free(img);
  return rval;
} 

void 
WriteData( const std::string fileName, 
  const double *data, const int nData) {

  std::cerr << "writing: " << fileName << std::endl;
  
  std::ofstream output( fileName.c_str() );
  
  for( int cData=0; cData<nData; cData++ ) {
    output << data[ cData ] << std::endl;
  }
  
  output.close();
    
}

void 
WriteData( const std::string fileName, 
  double **dataArray, const int nRows, const int nCols ) {

  std::cerr << "writing: " << fileName << std::endl;
  
  std::ofstream output( fileName.c_str() );
  
  for( int cRow=0; cRow<nRows; cRow++ ) {
    for( int cCol=0; cCol<nCols; cCol++ ) {
      output << dataArray[ cRow ][ cCol ] << "   ";
    }    
    output << std::endl;
  }
  
  output.close();
    
}

std::string 
GetFieldAndParameter( std::string field, std::string parameter ) {
  std::ostringstream output;
  output << field << ": \n" << "  " << parameter << std::endl << std::endl;
  return output.str();
}

std::string 
GetVersion() {
  return "1.0 Beta";
}

std::string 
GetCurrentDirectory() {
  
  const int nChars = 2000;
  char cwd[ nChars ];
  
  getcwd( cwd, nChars );
  
  return std::string( cwd );
}

std::string 
GetCommandArguments( int argc, char * argv[] ) {  

  std::ostringstream output;

  for( int cArg=0; cArg<argc; cArg++ ) {
    output << argv[ cArg ] << " ";
  }
  output << std::endl;
  
  return output.str();
}

std::string 
GetOperatingSystem() {  
  utsname uts;
  uname( &uts );    
  std::string name( uts.sysname );
  return name;
}

std::string
GetHostName() {
  utsname uts;
  uname( &uts );
  std::string name( uts.nodename );
  return name;
}

std::string
GetMachine() {
  utsname uts;
  uname( &uts );
  std::string machine( uts.machine );
  return machine;
}

std::string
GetEnvironmentVariable( std::string variable ) {
  char *result = std::getenv( variable.c_str() );
  std::string stringResult;
  if( result == NULL ) {
    stringResult = "not set";
  } else {
    stringResult = std::string( result );
  }
  return stringResult;
}

std::string 
GetFreeSurferHome() {
  return GetEnvironmentVariable( "FREESURFER_HOME" );
}

std::string 
GetSubjectsDirectory() {
  return GetEnvironmentVariable( "SUBJECTS_DIR" );
}

std::string
GetFileExtension( std::string fileName ) {

  // take in part from itk  
  std::string::size_type ext = fileName.rfind( '.' );
  std::string exts = fileName.substr( ext );
  if( exts == ".gz" ) {
    std::string::size_type dotpos = fileName.rfind( '.', ext-1 );
    if( dotpos != std::string::npos ) {
      exts = fileName.substr(dotpos);
    }
  }  
  
  return exts;
}

void PrintUsageError( std::string error ) {
  std::cerr << "\n*** usage error: " << error << std::endl << std::endl;
}

int main (int argc, char * argv[]) {

  PARSE_ARGS;
  
  // add the MGH/MGZ reader
  itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() ); 
    
  CommandUpdate::Pointer observer = CommandUpdate::New();
  
  if( argc <= 1 ) {
    std::cerr << "\nuse --help flag for usage\n" << std::endl;
    return EXIT_FAILURE;
  }

  if( !outputDirectory.empty() ) {
    observer->SetOutputDirectory( outputDirectory );
  } else {
    PrintUsageError( "must specify output directory" );
    return EXIT_FAILURE;
  }
  observer->SetLogFileName( "poistats.log" );

  observer->PostMessage( "-- Poistats " + GetVersion() + " --\n\n" );
  
  observer->PostMessage( GetFieldAndParameter( "Command", 
    GetCommandArguments( argc, argv ) ) );

  observer->PostMessage( GetFieldAndParameter( "FreeSurfer Home", 
    GetFreeSurferHome() ) );

  observer->PostMessage( GetFieldAndParameter( "Subjects Directory", 
    GetSubjectsDirectory() ) );

  observer->PostMessage( GetFieldAndParameter( "Current Working Directory", 
    GetCurrentDirectory() ) );
  
  observer->PostMessage( GetFieldAndParameter( "Operating System", 
    GetOperatingSystem() ) );

  observer->PostMessage( GetFieldAndParameter( "Host name", GetHostName() ) );

  observer->PostMessage( GetFieldAndParameter( "Machine", GetMachine() ) );
  
  std::string imageFileExtension = GetFileExtension( diffusionTensorImage );

  typedef TensorReaderStrategy::TensorPixelType TensorPixelType;
  typedef TensorReaderStrategy::TensorImageType TensorImageType;

  TensorReaderStrategy *tensorReader = NULL;
  if( !isSymmetricTensorData ){
    
    if( IsMGH( imageFileExtension ) ){
      tensorReader = new AsymmetricTensorVectorReaderStrategy();
    } else {
      tensorReader = new AsymmetricTensorReaderStrategy();
    }
    
  } else {
    tensorReader = new SymmetricTensorReaderStrategy();
  }
  tensorReader->SetObserver( observer );
  tensorReader->SetFileName( diffusionTensorImage );
  TensorImageType::Pointer tensors = tensorReader->GetTensors();
  
  delete tensorReader;
      
  typedef itk::Image< float, 3 > OutputImageType;
  typedef itk::PoistatsFilter< TensorImageType, OutputImageType > 
    PoistatsFilterType;
  PoistatsFilterType::Pointer poistatsFilter = PoistatsFilterType::New();
  
  poistatsFilter->AddObserver( itk::AnyEvent(), observer );

  // this is a hack to get the directions loaded corrently.  There's a bug
  // in the itk nifti reader that sets the direction cosines to the wrong
  // to the canonical directions rather than the true calculated directions
  if( IsNifti( imageFileExtension ) ) {

    mat44 transform = GetNiftiTransform( diffusionTensorImage.c_str() );
  
    TensorImageType::DirectionType direction = tensors->GetDirection();
    for( int cRow=0; cRow<3; cRow++ ) {
      for( int cCol=0; cCol<3; cCol++ ) {
        direction( cRow, cCol ) = transform.m[ cRow ][ cCol ] / 
          tensors->GetSpacing()[ cRow ];
      }
    }
        
    tensors->SetDirection( direction );
    
  }
  
  std::cerr << "direction: \n" << tensors->GetDirection();
  poistatsFilter->SetInput( tensors );

  observer->PostMessage( "reading seed volume...\n" );
  // read seed volume
  typedef itk::ImageFileReader< PoistatsFilterType::SeedVolumeType > SeedReaderType;
  SeedReaderType::Pointer seedReader = SeedReaderType::New();
  seedReader->SetFileName( seedRegions );
  try { 
    seedReader->Update();
  } catch( itk::ExceptionObject & excp ) {
    std::ostringstream output;
    output << "Error reading the series." << std::endl << excp << std::endl;
    observer->PostErrorMessage( output.str() );
    return EXIT_FAILURE;
  }
  poistatsFilter->SetSeedVolume( seedReader->GetOutput() );
  
  // set seed for random number generator
  std::srand( ( unsigned ) time( 0 ) ); 
  poistatsFilter->SetRandomSeed( std::rand() );

  // set polarity
  PoistatsFilterType::MatrixType polarity( 3, 3 );
  polarity[ 0 ][ 0 ] = 1;
  polarity[ 0 ][ 1 ] = 1;
  polarity[ 0 ][ 2 ] = -1;

  polarity[ 1 ][ 0 ] = 1;
  polarity[ 1 ][ 1 ] = 1;
  polarity[ 1 ][ 2 ] = -1;

  polarity[ 2 ][ 0 ] = -1;
  polarity[ 2 ][ 1 ] = -1;
  polarity[ 2 ][ 2 ] = 1;
  poistatsFilter->SetPolarity( polarity );

  poistatsFilter->SetNumberOfControlPoints( numberOfControlPoints );
  poistatsFilter->SetInitialSigma( sigma );

  observer->PostMessage( "reading sampling volume...\n" );

  // read sampling volume
  typedef itk::ImageFileReader< PoistatsFilterType::SamplingVolumeType > SamplingReaderType;
  SamplingReaderType::Pointer samplingReader = SamplingReaderType::New();
  samplingReader->SetFileName( sample );
  try { 
    samplingReader->Update();
  } catch( itk::ExceptionObject & excp ) {
    std::ostringstream output;
    output << "Error reading the series." << std::endl << excp << std::endl;
    observer->PostErrorMessage( output.str() );
    return EXIT_FAILURE;
  }
  poistatsFilter->SetSamplingVolume( samplingReader->GetOutput() );
  
  // set number of sample points
  poistatsFilter->SetNumberOfSamplePoints( numberOfSamplePoints );

  // read mask volume if it exists
  if( mask.size() != 0 ) {
    observer->PostMessage( "reading mask...\n" );
    typedef itk::ImageFileReader< PoistatsFilterType::MaskVolumeType > MaskReaderType;
    MaskReaderType::Pointer maskReader = MaskReaderType::New();
    maskReader->SetFileName( mask );
    try { 
      maskReader->Update();
    } catch( itk::ExceptionObject & excp ) {
      std::ostringstream output;
      output << "Error reading the series." << std::endl << excp << std::endl;
      observer->PostErrorMessage( output.str() );
      return EXIT_FAILURE;
    }
    poistatsFilter->SetMaskVolume( maskReader->GetOutput() );
  }

  // seeds the seeds to be used if they exist  
  for( std::vector< int >::iterator seedIt = seedRegionsToUse.begin();
    seedIt != seedRegionsToUse.end(); seedIt++ ) {
  
    int seedValueToUse = ( *seedIt );
    poistatsFilter->SetNextSeedValueToUse( seedValueToUse );      
  }
  
  if( exchangeProbability > 0 ) {
    poistatsFilter->SetReplicaExchangeProbability( exchangeProbability );
  }

  if( timeConstant > 0 ) {
    poistatsFilter->SetSigmaTimeConstant( timeConstant );
  }
  
  if( gamma > 0 ) {
    poistatsFilter->SetPointsToImageGamma( gamma );
  }

  if( numberOfReplicas > 0 ) {
    poistatsFilter->SetNumberOfReplicas( numberOfReplicas );
  }

  // compute the poi
  try { 
    poistatsFilter->Update();
  } catch( itk::ExceptionObject & excp ) {
    std::ostringstream output;
    output << "Error thrown in poistats filter." << std::endl << excp << std::endl;
    observer->PostErrorMessage( output.str() );
    return EXIT_FAILURE;
  }
 
 // this overrides the default output file format
  if( isOutputNii ){
    imageFileExtension = ".nii";
  }
  
  typedef itk::ImageFileWriter< PoistatsFilterType::OutputImageType > WriterType;  
  WriterType::Pointer writer = WriterType::New();

  // write aggregate densities
  std::string densityFileName = (std::string)outputDirectory + 
    (std::string)"/PathDensity" + imageFileExtension;
  OutputImageType::Pointer pathDensity = 
    poistatsFilter->GetOutput( PoistatsFilterType::PATH_DENSITY_OUTPUT );
  writer->SetInput( pathDensity );
  writer->SetFileName( densityFileName.c_str() );
  
  observer->PostMessage( "writing: " + densityFileName + "\n" );  
  writer->Update();  

  std::string optimalDensityFileName = (std::string)outputDirectory + 
    (std::string)"/OptimalPathDensity" + imageFileExtension;
  OutputImageType::Pointer optimalPathDensity = 
    poistatsFilter->GetOutput( PoistatsFilterType::OPTIMAL_PATH_DENSITY_OUTPUT );
  writer->SetInput( optimalPathDensity );
  writer->SetFileName( optimalDensityFileName.c_str() );

  observer->PostMessage( "writing: " + optimalDensityFileName + "\n" );  
  writer->Update();
  
  PoistatsFilterType::MatrixType finalPath = poistatsFilter->GetFinalPath();
  const std::string finalPathFileName( (std::string)outputDirectory + 
    (std::string)"/OptimalPath.txt" );
  WriteData( finalPathFileName, finalPath.data_array(), 
    finalPath.rows(), finalPath.cols() );

  PoistatsFilterType::ArrayType fa = poistatsFilter->GetSamples();  

  const std::string faFileName( (std::string)outputDirectory + 
    (std::string)"/OptimalPathSamples.txt" );
  WriteData( faFileName, fa.data_block(), fa.size() );
  
  PoistatsFilterType::ArrayType finalPathProbabilities = 
    poistatsFilter->GetFinalPathProbabilities();
  const std::string pathProbabilitiesFileName( (std::string)outputDirectory + 
    (std::string)"/OptimalPathProbabilities.txt" );
  WriteData( pathProbabilitiesFileName, 
    finalPathProbabilities.data_block(), finalPathProbabilities.size() );
    
  return EXIT_SUCCESS;

}
