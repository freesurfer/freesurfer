#include <iostream>
#include <string>

extern "C"
{
#include "matrix.h"
}

#include "ui/CommandParser.h"
#include "ui/FreeSurferExecutable.h"

#include "datamodel/utils/InitializePath.h"

/** This is needed by the freesurfer utils library */
char *Progname;

/**
 * C++ replacement of the matlab Poistats.
 */
class InitializePathExe : public FreeSurferExecutable {

  public:
    
    // string indicating a flag
    static const std::string FLAG_DELIMITER;
    
    // principle eigenvector input
    static const std::string FLAG_EIGENVECTOR;
    
    // seed volume input
    static const std::string FLAG_SEEDS;

    // seed values to connect    
    static const std::string FLAG_SEED_VALUES;

    // output directory
    static const std::string FLAG_OUTPUT_DIRECTORY;

    InitializePathExe( int inArgs, char ** iaArgs );
    ~InitializePathExe();

    /**
     * Fills in the arguments and returns true of the arguments can be filled.
     */
    bool FillArguments();    
    
    void Run();
    
    MRI* ReadMRI( const char* fileName );
  
  private:

    const char *m_EigenVectorFileName;
    const char *m_SeedVolumeFileName;
    std::vector< int > m_SeedValues;   
    const char *m_OutputDir;

};

const std::string InitializePathExe::FLAG_DELIMITER = "--";

const std::string InitializePathExe::FLAG_EIGENVECTOR = 
  InitializePathExe::FLAG_DELIMITER + "ev";

const std::string InitializePathExe::FLAG_SEEDS = 
  InitializePathExe::FLAG_DELIMITER + "seeds";

const std::string InitializePathExe::FLAG_SEED_VALUES = 
  InitializePathExe::FLAG_DELIMITER + "seednums";

const std::string InitializePathExe::FLAG_OUTPUT_DIRECTORY = 
  InitializePathExe::FLAG_DELIMITER + "o";

InitializePathExe::InitializePathExe( int inArgs, char ** iaArgs ) : 
  FreeSurferExecutable( inArgs, iaArgs ) {

  m_EigenVectorFileName = NULL;
  m_SeedVolumeFileName = NULL;
  m_OutputDir = NULL;

  SetName( "dmri_initialize_path", "find starting path for input to poistats" );  

  SetNextRequiredArgument( FLAG_EIGENVECTOR, "eigvec", 
    "Principle eigenvector input", 
    "eigvec1.nii", "must specify an eigenvector volume" );

  SetNextRequiredArgument( FLAG_SEEDS, "seeds", 
    "Seed volume", 
    "seeds.nii", "must specify a seed volume" );

  SetNextRequiredArgument( FLAG_OUTPUT_DIRECTORY, "output", 
    "Output directory", 
    "output/", "must specify an output directory" );

  SetNextOptionalArgument( FLAG_SEED_VALUES, "seednumvalue", 
    "Use <seednumvalue> to define seed region.  If nothing is specified, then all seeds are used.  More than one must be available.  Eg, --seednums 1,2" );


  std::string output = "";
  output = output + 
    "InitialPath - volume of initial points" + "\n   " +
    "InitialPath.txt - coordinates of initial path";

  SetOutput( output );
  
  SetVersion( "0.1" );
  SetBugEmail( "martinos-tech@yahoogroups.com" );
    
}

InitializePathExe::~InitializePathExe() {  

  if( m_EigenVectorFileName != NULL ) {
    delete m_EigenVectorFileName; 
  }
  
  if( m_SeedVolumeFileName != NULL ) {
    delete m_SeedVolumeFileName; 
  }
  
  if( m_OutputDir != NULL ) {
    delete m_OutputDir; 
  }
  
}

bool
InitializePathExe::FillArguments() {
  bool isFilled = false;

  try {
    std::string *requiredArguments = GetRequiredArguments();  

    m_EigenVectorFileName = requiredArguments[0].c_str();
    m_SeedVolumeFileName = requiredArguments[1].c_str();
    m_OutputDir = requiredArguments[2].c_str();
    
    isFilled = true;    
  } catch(...) {
    isFilled = false;
  }

  if( isFilled ) {  
    m_SeedValues = m_Parser->GetArgumentIntVector( FLAG_SEED_VALUES.c_str() );
  }
  
  return isFilled;
}

void
InitializePathExe::Run() {
  
  InitializePath initializePath;
  
  // read eigen vectors
  MRI* eigenVectors = this->ReadMRI( m_EigenVectorFileName );
  initializePath.SetEigenVectors( eigenVectors );
  
  // read seed volume
  MRI* seedVolume = this->ReadMRI( m_SeedVolumeFileName );
  initializePath.SetSeedVolume( seedVolume );
  
  // set seed values to connect
  initializePath.SetSeedValues( &m_SeedValues );
  
  // now that all the inputs are set, calculate the initial path
  initializePath.CalculateInitialPath();
  
  // get the resulting pathway
  MATRIX *path = initializePath.GetInitialPath();
  
// TODO: removing this temporarily  
//  // save the pathway to the output directory
//  const std::string pathFileName( (std::string)m_OutputDir + 
//    (std::string)"/InitialPath.txt" );
//  WriteData( pathFileName, path->data, path->rows, path->cols );
  
  if( path != NULL ) {
    MatrixFree( &path );
    path = NULL;
  }

  MRIfree( &eigenVectors );
  MRIfree( &seedVolume );
}

MRI* 
InitializePathExe::ReadMRI( const char* fileName ) {

  // read the eigenvectors
  char* nonConstFileName = strdup( fileName );
  MRI* volume = ::MRIread( nonConstFileName );
  free( nonConstFileName );

  if ( volume == NULL ) {
    // TODO: throw an error
    std::cerr << "mri read failed..." << std::endl;
  }
  
  return volume;
}

int main( int argc, char ** argv ) {
  
  // this is needed by the freesurfer utils library
  Progname = argv[0];
  
  FreeSurferExecutable *exe = new InitializePathExe( argc, argv );
  
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
  
  delete exe;
  exe = NULL;
  
  return 0;
}
