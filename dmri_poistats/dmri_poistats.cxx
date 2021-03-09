/**
 * @brief Finds the most probable path between two seed regions.
 *
 */
/*
 * Original Author: Dennis Jen, Dave Tuch (matlab prototype) 
 *    $Author $
 *    $Date $
 *    $Revision $
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.eduf
 *
 */

#include <iostream>
#include <string>
#include <cstdlib>
#include <cstddef>

// for reading in the cosine directions because itk doesn't do it properly
#include <nifti1_io.h>

#include <itkImageSeriesReader.h>
#include <itkImageFileWriter.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkNiftiImageIO.h>
#include <itkOrientedImage.h>

#include "itkMGHImageIOFactory.h" // itkio
#include "io/SymmetricTensorReaderStrategy.h"
#include "io/AsymmetricTensorReaderStrategy.h"
#include "io/AsymmetricTensorVectorReaderStrategy.h"

#include "datamodel/utils/itkPoistatsFilter.h"

#include "datamodel/events/CommandUpdate.h"

#include "ui/CommandParser.h" //itkutils
#include "ui/FreeSurferExecutable.h" //itkutils

/** This is needed by the freesurfer utils library */
char *Progname;

/**
 * C++ replacement of the matlab Poistats.
 */
class Poistats : public FreeSurferExecutable {

  public:
  
    static const std::string FLAG_DELIMITER;
  
    // required
    static const std::string FLAG_INPUT_STEM;
    static const std::string FLAG_OUTPUT_DIRECTORY;
    static const std::string FLAG_SEEDS;
    static const std::string FLAG_SAMPLE_STEM;
    static const std::string FLAG_NUM_CONTROL_POINTS;
    static const std::string FLAG_SIGMA;
  
    // optional
    static const std::string FLAG_MASK_STEM;
    static const std::string FLAG_NUM_SAMPLES;
    static const std::string FLAG_SEED_VALUES;    
    static const std::string FLAG_NUM_REPLICAS;      
    static const int DEFAULT_NUM_SEEDS = 100;    
    static const std::string FLAG_REPLICA_EXCHANGE_PROBABILITY;
    static const std::string FLAG_SIGMA_TIME_CONSTANT;    
    static const std::string FLAG_POINTS_TO_IMAGE_GAMMA;
    static const std::string FLAG_IS_SYMMETRIC_DATA;
    static const std::string FLAG_IS_OUTPUT_NII;
    static const std::string FLAG_EIGENVECTOR_STEM;
    static const std::string FLAG_FIELD_LINE;
    static const std::string FLAG_SHOULD_WRITE_BEST_PATHS;

    Poistats( int inArgs, char ** iaArgs );
    ~Poistats();

    /**
     * Fills in the arguments and returns true of the arguments can be filled.
     */
    bool FillArguments();    
    void Run();
  
  private:

    // required arguments
    const char *m_InputStem;
    const char *m_OutputDir;
    const char *m_Seeds;
    const char *m_SampleStem;  
    int m_nControlPoints;
    int m_Sigma;
    
    // optional arguments
    char *m_MaskStem;
    int m_nSamples;
    std::vector< int > m_SeedValues;
    int m_nReplicas;
  
    double m_ReplicaExchangeProbability;  
    double m_SigmaTimeConstant;      
    double m_PointsToImageGamma;
  
    bool m_IsSymmetricData;
    bool m_IsOutputNii;
    
    char *m_EigenVectorStem;
    
    double m_FieldLineRadius;
    
    bool m_ShouldWriteBestPaths;

    bool IsNifti( std::string fileExtension );
    bool IsMGH( std::string fileExtension );
    
    mat44 GetNiftiTransform( const char *filename );

};

const std::string Poistats::FLAG_DELIMITER = "--";

// required flags
const std::string Poistats::FLAG_INPUT_STEM = Poistats::FLAG_DELIMITER + "i";
const std::string Poistats::FLAG_OUTPUT_DIRECTORY = 
  Poistats::FLAG_DELIMITER + "o";
const std::string Poistats::FLAG_SEEDS = 
  Poistats::FLAG_DELIMITER + "seeds";
const std::string Poistats::FLAG_SAMPLE_STEM = 
  Poistats::FLAG_DELIMITER + "sample";
const std::string Poistats::FLAG_NUM_CONTROL_POINTS = 
  Poistats::FLAG_DELIMITER + "nc";
const std::string Poistats::FLAG_SIGMA = Poistats::FLAG_DELIMITER + "sigma";

// optional flags
const std::string Poistats::FLAG_MASK_STEM = Poistats::FLAG_DELIMITER + "m";
const std::string Poistats::FLAG_NUM_SAMPLES = Poistats::FLAG_DELIMITER + "ns";
const std::string Poistats::FLAG_SEED_VALUES = 
  Poistats::FLAG_DELIMITER + "seednums";
const std::string Poistats::FLAG_NUM_REPLICAS = 
  Poistats::FLAG_DELIMITER + "nreplicas";
const std::string Poistats::FLAG_IS_SYMMETRIC_DATA = 
  Poistats::FLAG_DELIMITER + "symmetric";
const std::string Poistats::FLAG_IS_OUTPUT_NII = 
  Poistats::FLAG_DELIMITER + "nii";
const std::string Poistats::FLAG_EIGENVECTOR_STEM = 
  Poistats::FLAG_DELIMITER + "eigvec";
const std::string Poistats::FLAG_FIELD_LINE = 
  Poistats::FLAG_DELIMITER + "fieldline";
const std::string Poistats::FLAG_SHOULD_WRITE_BEST_PATHS =
  Poistats::FLAG_DELIMITER + "writepaths";

// these are new additions for playing with the parameter space
const std::string Poistats::FLAG_REPLICA_EXCHANGE_PROBABILITY = 
  Poistats::FLAG_DELIMITER + "exchangeprob";
const std::string Poistats::FLAG_SIGMA_TIME_CONSTANT = 
  Poistats::FLAG_DELIMITER + "timeconst";
const std::string Poistats::FLAG_POINTS_TO_IMAGE_GAMMA = 
  Poistats::FLAG_DELIMITER + "gamma";

Poistats::Poistats( int inArgs, char ** iaArgs ) : 
  FreeSurferExecutable( inArgs, iaArgs ) {
  SetName( "dmri_poistats", "find optimal path in tensor volume (**beta version with limited support)" );  

  SetNextRequiredArgument( FLAG_INPUT_STEM, "dtensorinstem", 
    "Diffusion tensor input", 
    "dtensor.nii", "must specify a dtensor input filename" );
  SetNextRequiredArgument( FLAG_OUTPUT_DIRECTORY, "outdir", 
    "Output directory.",
    "poistats", "must specify an output directory" );
  SetNextRequiredArgument( FLAG_SEEDS, "seedstem", 
    "Volume containing numerical labels to use as seed regions.", 
    "seedvol.nii", "must specify a seed volume" );
  SetNextRequiredArgument( FLAG_SAMPLE_STEM, "samplestem", 
    "Instem for volume to sample. For example: fa, trace", "fa.nii", 
    "must specify a sampling volume" );
  SetNextRequiredArgument( FLAG_NUM_CONTROL_POINTS, "ncontrolpoints", 
    "Number of control points used to describe path. Number should be approximately the number of 'turns' in the path. Almost always 1 or 2.", 
    "2", "must specify number of control points" );
  SetNextRequiredArgument( FLAG_SIGMA, "sigmasize", 
    "Search distance for path control points in units voxels. The search distance should be approximately the distance between the midpoint of the seed path and the target location", 
    "10", "must specify search distance" );

  SetNextOptionalArgument( FLAG_MASK_STEM, "maskstem", 
    "Instem for mask. The path will not be allowed to contact the mask. The mask can be used to specify invalid regions, e.g., CSF" );
  SetNextOptionalArgument( FLAG_NUM_SAMPLES, "nsamplepoints", 
    "Number of points to sample along path from sample volume. For example, --ns 100 will sample 100 values along the path. Default: 100" );
  SetNextOptionalArgument( FLAG_SEED_VALUES, "seednumvalue", 
    "Use <seednumvalue> to define seed region. Eg, --seednums 1,2" );
  SetNextOptionalArgument( FLAG_NUM_REPLICAS, "nreplicas", 
    "Use <nreplicas> to specify the number of replicas.  For example, --nreplicas 100 will spawn 100 replicas. Default: 100" );

  SetNextOptionalArgument( FLAG_REPLICA_EXCHANGE_PROBABILITY, "exchangeprob", 
    "Replica exchange probability.  Default: 0.05" );
  SetNextOptionalArgument( FLAG_SIGMA_TIME_CONSTANT, "timeconst", 
    "Sigma time constant.  Default: 200." );
  SetNextOptionalArgument( FLAG_POINTS_TO_IMAGE_GAMMA, "gamma", 
    "Points to image gamma.  Default: 0.5." );
    
  SetNextOptionalArgument( FLAG_IS_SYMMETRIC_DATA, "symmetric", 
    "Is the tensor input stored symmetric (6 rather than 9 component storage)?  Default: false (not stored symmetrically)." );
  SetNextOptionalArgument( FLAG_IS_OUTPUT_NII, "nii", 
    "Should the output be stored as nifti?  Default: false" );
  SetNextOptionalArgument( FLAG_EIGENVECTOR_STEM, "eigvec", 
    "Should the paths be initialized by following the principle eigenvector?  Specify the eigenvector volume following this flag." );
  SetNextOptionalArgument( FLAG_FIELD_LINE, "fieldline", 
    "Should the paths be initialized by field lines?  Specify the maximum amplitude." );
  SetNextOptionalArgument( FLAG_SHOULD_WRITE_BEST_PATHS, "writepaths",
    "TESTING -- Should the best paths from each replicas be written out?  This will save to InitialPaths.txt for the time being" );

  std::string output = "";
  output = output + 
    "PathDensity.nii - path probability density" + "\n   " +
    "OptimalPathDensity.nii - probability density of optimal path" + "\n\n   " +

    "OptimalPathSamples.txt - values of sample volume along optimal path" + "\n   " +
    "OptimalPathProbabilities.txt - probability values along optimal path" + "\n   " +
    "OptimalPath.txt - coordinates of optimal path";

  SetOutput( output );
  
  SetVersion( "1.0 beta" );
  SetBugEmail( "freesurfer@nmr.mgh.harvard.edu" );
}

Poistats::~Poistats() {  
}

bool
Poistats::FillArguments() {
  bool isFilled = false;

  try {
    std::string *requiredArguments = GetRequiredArguments();  
    m_InputStem = requiredArguments[0].c_str();
    m_OutputDir = requiredArguments[1].c_str();
    m_Seeds = requiredArguments[2].c_str();
    m_SampleStem = requiredArguments[3].c_str();
    
    m_nControlPoints = atoi( requiredArguments[4].c_str() );
    m_Sigma = atoi( requiredArguments[5].c_str() );
    
    isFilled = true;    
  } catch(...) {
    isFilled = false;
  }    
  
  if( isFilled ) {
    // optional parameters
    m_MaskStem = m_Parser->GetArgument( FLAG_MASK_STEM.c_str() );
  
    m_nSamples = m_Parser->GetArgumentInt( FLAG_NUM_SAMPLES.c_str() );
    if( m_nSamples <= 0) {
      m_nSamples = DEFAULT_NUM_SEEDS;
      
      // TODO: use the observer
      std::cout << "INFO: setting nsamplepoints=" << m_nSamples << std::endl;
    }
    
    m_SeedValues = m_Parser->GetArgumentIntVector( FLAG_SEED_VALUES.c_str() );

    m_nReplicas = m_Parser->GetArgumentInt( FLAG_NUM_REPLICAS.c_str() );
    
    m_ReplicaExchangeProbability = m_Parser->GetArgumentDouble( 
      FLAG_REPLICA_EXCHANGE_PROBABILITY.c_str() );

    m_SigmaTimeConstant = m_Parser->GetArgumentDouble( 
      FLAG_SIGMA_TIME_CONSTANT.c_str() );

    m_PointsToImageGamma = m_Parser->GetArgumentDouble( 
      FLAG_POINTS_TO_IMAGE_GAMMA.c_str() );

    m_IsSymmetricData = m_Parser->GetArgumentBoolean( 
      FLAG_IS_SYMMETRIC_DATA.c_str() );

    m_IsOutputNii = m_Parser->GetArgumentBoolean( FLAG_IS_OUTPUT_NII.c_str() );

    m_EigenVectorStem = m_Parser->GetArgument( FLAG_EIGENVECTOR_STEM.c_str() );

    m_FieldLineRadius = m_Parser->GetArgumentDouble( FLAG_FIELD_LINE.c_str() );
    
    m_ShouldWriteBestPaths = m_Parser->GetArgumentBoolean( FLAG_SHOULD_WRITE_BEST_PATHS.c_str() );
    
  }
      
  // TODO: delete requiredArguments  
  
  return isFilled;
}

bool
Poistats::IsNifti( std::string fileExtension ) {

  bool isNifti = false;
  
  if( fileExtension == ".nii" || fileExtension == ".nii.gz" ) {
    isNifti = true;
  }
  
  return isNifti;
}

bool
Poistats::IsMGH( std::string fileExtension ) {

  bool isMGH = false;
  
  if( fileExtension == ".mgh" || fileExtension == ".mgz" ) {
    isMGH = true;
  }
  
  return isMGH;
}

mat44
Poistats::GetNiftiTransform( const char *filename ) {
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
Poistats::Run() {
  // add the MGH/MGZ reader
  itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() ); 
    
  CommandUpdate::Pointer observer = CommandUpdate::New();
    
  if( m_OutputDir != NULL ) {
    observer->SetOutputDirectory( m_OutputDir );

// TODO: create output directory
//      err = mkdir(outdir,0777);
//      if (err != 0 && errno != EEXIST) {
//        printf("ERROR: creating directory %s\n",outdir);
//        perror(NULL);
//        exit(1);
//      } 

  } else {
    PrintUsageError( "must specify output directory" );
    return;
  }
  observer->SetLogFileName( "poistats.log" );

  observer->PostMessage( "-- Poistats " + GetVersion() + " --\n\n" );
  
  observer->PostMessage( GetFieldAndParameter( "Command", 
    m_CommandLineArguments ) );

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
  
  std::string imageFileExtension = GetFileExtension( m_InputStem );

  typedef TensorReaderStrategy::TensorPixelType TensorPixelType;
  typedef TensorReaderStrategy::TensorImageType TensorImageType;

  TensorReaderStrategy *tensorReader = NULL;
  if( !m_IsSymmetricData ){
    
    if( IsMGH( imageFileExtension ) ){
      tensorReader = new AsymmetricTensorVectorReaderStrategy();
    } else {
      tensorReader = new AsymmetricTensorReaderStrategy();
    }
    
  } else {
    tensorReader = new SymmetricTensorReaderStrategy();
  }
  tensorReader->SetObserver( observer );
  tensorReader->SetFileName( m_InputStem );
  TensorImageType::Pointer tensors = tensorReader->GetTensors();
  
  delete tensorReader;
      
  typedef itk::OrientedImage< float, 3 > OutputImageType;
  typedef itk::PoistatsFilter< TensorImageType, OutputImageType > 
    PoistatsFilterType;
  PoistatsFilterType::Pointer poistatsFilter = PoistatsFilterType::New();
  
  poistatsFilter->AddObserver( itk::AnyEvent(), observer );

  // this is a hack to get the directions loaded corrently.  There's a bug
  // in the itk nifti reader that sets the direction cosines to the wrong
  // to the canonical directions rather than the true calculated directions
  if( IsNifti( imageFileExtension ) ) {

    mat44 transform = GetNiftiTransform( m_InputStem );
  
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
  seedReader->SetFileName( m_Seeds );
  try { 
    seedReader->Update();
  } catch( itk::ExceptionObject & excp ) {
    std::ostringstream output;
    output << "Error reading the series." << std::endl << excp << std::endl;
    observer->PostErrorMessage( output.str() );
    return;
  }
  poistatsFilter->SetSeedVolume( seedReader->GetOutput() );
  
  // read in the seeds using the MGH reader
  MRI *mghSeedVolume = FreeSurferExecutable::ReadMRI( m_Seeds );
  
  // set the seeds
  poistatsFilter->SetMghSeeds( mghSeedVolume );

  // set seed for random number generator
  std::srand( ( unsigned ) time( 0 ) ); 
  poistatsFilter->SetRandomSeed( std::rand() );

  // set polarity -- specific for the scanner that was used to collect the data
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

  poistatsFilter->SetNumberOfControlPoints( m_nControlPoints );
  poistatsFilter->SetInitialSigma( m_Sigma );

  observer->PostMessage( "reading sampling volume...\n" );

  // read sampling volume
  typedef itk::ImageFileReader< PoistatsFilterType::SamplingVolumeType > SamplingReaderType;
  SamplingReaderType::Pointer samplingReader = SamplingReaderType::New();
  samplingReader->SetFileName( m_SampleStem );
  try { 
    samplingReader->Update();
  } catch( itk::ExceptionObject & excp ) {
    std::ostringstream output;
    output << "Error reading the series." << std::endl << excp << std::endl;
    observer->PostErrorMessage( output.str() );
    return;
  }
  poistatsFilter->SetSamplingVolume( samplingReader->GetOutput() );
  
  // set number of sample points
  poistatsFilter->SetNumberOfSamplePoints( m_nSamples );

  // read mask volume if it exists
  if( m_MaskStem != NULL ) {
    observer->PostMessage( "reading mask...\n" );
    typedef itk::ImageFileReader< PoistatsFilterType::MaskVolumeType > MaskReaderType;
    MaskReaderType::Pointer maskReader = MaskReaderType::New();
    maskReader->SetFileName( m_MaskStem );
    try { 
      maskReader->Update();
    } catch( itk::ExceptionObject & excp ) {
      std::ostringstream output;
      output << "Error reading the series." << std::endl << excp << std::endl;
      observer->PostErrorMessage( output.str() );
      return;
    }
    poistatsFilter->SetMaskVolume( maskReader->GetOutput() );
  }

  // seeds the seeds to be used if they exist  
  if( m_SeedValues.size() == 1 ) {
    std::ostringstream output;
    output << "If specifying seed number, you must specify at least two labels." << std::endl;
    observer->PostErrorMessage( output.str() );
    return;
  } else {
    for( std::vector< int >::iterator seedIt = m_SeedValues.begin();
      seedIt != m_SeedValues.end(); seedIt++ ) {        
      int seedValueToUse = ( *seedIt );
      poistatsFilter->SetNextSeedValueToUse( seedValueToUse );      
    }
  } 
  
  if( m_ReplicaExchangeProbability > 0 ) {
    poistatsFilter->SetReplicaExchangeProbability( 
      m_ReplicaExchangeProbability );
  }

  if( m_SigmaTimeConstant > 0 ) {
    poistatsFilter->SetSigmaTimeConstant( m_SigmaTimeConstant );
  }
  
  if( m_PointsToImageGamma > 0 ) {
    poistatsFilter->SetPointsToImageGamma( m_PointsToImageGamma );
  }

  if( m_nReplicas > 0 ) {
    poistatsFilter->SetNumberOfReplicas( m_nReplicas );
  }
  
  MRI *mghEigenVectors = NULL;
  
  // give the field line initialization priority over the eigenvector
  // initialization  
  if( m_FieldLineRadius > 0 ) {
    observer->PostMessage( "field line initialization...\n" );    
    poistatsFilter->SetFieldLineRadius( m_FieldLineRadius );
    poistatsFilter->SetUsingFieldLineInitialization();
  } else if( m_EigenVectorStem != NULL ) {
    observer->PostMessage( "reading eigenvectors...\n" );
    mghEigenVectors = FreeSurferExecutable::ReadMRI( m_EigenVectorStem );
    poistatsFilter->SetMghEigenVectors( mghEigenVectors );
    poistatsFilter->SetUsingEigenVectorInitialization();
  } else {
    poistatsFilter->SetUsingNormalInitialization();
  }

  // compute the poi
  try { 
    poistatsFilter->Update();
  } catch( itk::ExceptionObject & excp ) {
    std::ostringstream output;
    output << "Error thrown in poistats filter." << std::endl << excp << std::endl;
    observer->PostErrorMessage( output.str() );
    return;
  }
 
 // this overrides the default output file format
  if( m_IsOutputNii ){
    imageFileExtension = ".nii";
  } else {
    // for some reason, the nii writer isn't writing the directions out 
    // correctly, so by default, we're going to write out as mgz for the time 
    // being
    imageFileExtension = ".mgz";
  }
  
  typedef itk::ImageFileWriter< PoistatsFilterType::OutputImageType > WriterType;  
  WriterType::Pointer writer = WriterType::New();

  // write aggregate densities
  std::string densityFileName = (std::string)m_OutputDir + 
    (std::string)"/PathDensity" + imageFileExtension;
  OutputImageType::Pointer pathDensity = 
    poistatsFilter->GetOutput( PoistatsFilterType::PATH_DENSITY_OUTPUT );
  writer->SetInput( pathDensity );
  writer->SetFileName( densityFileName.c_str() );
  
  observer->PostMessage( "writing: " + densityFileName + "\n" );  
  writer->Update();  

  std::string optimalDensityFileName = (std::string)m_OutputDir + 
    (std::string)"/OptimalPathDensity" + imageFileExtension;
  OutputImageType::Pointer optimalPathDensity = 
    poistatsFilter->GetOutput( PoistatsFilterType::OPTIMAL_PATH_DENSITY_OUTPUT );
  writer->SetInput( optimalPathDensity );
  writer->SetFileName( optimalDensityFileName.c_str() );

  observer->PostMessage( "writing: " + optimalDensityFileName + "\n" );  
  writer->Update();
  
  PoistatsFilterType::MatrixType finalPath = poistatsFilter->GetFinalPath();
  const std::string finalPathFileName( (std::string)m_OutputDir + 
    (std::string)"/OptimalPath.txt" );
  WriteData( finalPathFileName, finalPath.data_array(), 
    finalPath.rows(), finalPath.cols() );

  PoistatsFilterType::ArrayType fa = poistatsFilter->GetSamples();  

  const std::string faFileName( (std::string)m_OutputDir + 
    (std::string)"/OptimalPathSamples.txt" );
  WriteData( faFileName, fa.data_block(), fa.size() );
  
  PoistatsFilterType::ArrayType finalPathProbabilities = 
    poistatsFilter->GetFinalPathProbabilities();
  const std::string pathProbabilitiesFileName( (std::string)m_OutputDir + 
    (std::string)"/OptimalPathProbabilities.txt" );
  WriteData( pathProbabilitiesFileName, 
    finalPathProbabilities.data_block(), finalPathProbabilities.size() );
    
  // if we should write out the best paths, then get them and write them
  if( m_ShouldWriteBestPaths ) {
    // save the pathway to the output directory
    const std::string pathFileName( (std::string)m_OutputDir + 
      (std::string)"/InitialPath.txt" );
      
    PoistatsFilterType::MatrixListType paths = poistatsFilter->GetBestTrialPaths();
    
    bool isFirst = true;
    
    // iterate through the paths and write append them
    for( PoistatsFilterType::MatrixListType::iterator replicaIterator = paths.begin(); 
      replicaIterator != paths.end(); replicaIterator++ ) {

      PoistatsFilterType::MatrixPointer path = *replicaIterator;
            
      if( isFirst ) {

        // overwrite previous enteries
        WriteData( pathFileName, path->data_array(), path->rows(), path->cols() );
        isFirst = false;

      } else {
        
        // only write out the delimiter if this is not the first path
        const float delimiter = -999;
        float delimiterPoint[] = { delimiter, delimiter, delimiter };
        WriteDataAppend( pathFileName, delimiterPoint, 1, 3 );

        // write the pathway appended
        WriteDataAppend( pathFileName, path->data_array(), path->rows(), path->cols() );      

      }
      
    }
    
  }
    
  if( mghEigenVectors != NULL ) {
    MRIfree( &mghEigenVectors );
  }
  
  if( mghSeedVolume != NULL ) {
    MRIfree( &mghSeedVolume );
  }
  
    
  return;
}

int main( int argc, char ** argv ) {

  // this is needed by the freesurfer utils library
  Progname = argv[0];
  
  FreeSurferExecutable *exe = new Poistats( argc, argv );
  
  bool shouldPrintHelp = false;

  if( argc == 1 ) {
    shouldPrintHelp = true;
  } else if( argc == 2 ) {
    if( argv[ 1 ] == std::string( "--help" ) ) {
      shouldPrintHelp = true;
    }
  }
  
  if( shouldPrintHelp ) {
    exe->PrintHelp();
  } else {
    bool isFilled = exe->FillArguments();
    if( isFilled ) {
      exe->Run();
    } else {
      std::cerr << "\nuse --help for usage" << std::endl;
    }
  }
  
  return 0;
}

