// this needs to be included or a vxl linking error emerges
#include <itkBSplineInterpolateImageFunction.h>

#include <iostream>
#include <string>



#include "matrix.h"


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
    
    // number of times to repeat initialization
    static const std::string FLAG_NUMBER_OF_REPLICAS;
    
    // maximum distance 
    static const std::string FLAG_RADIUS;

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
    int m_NumberOfReplicas;
    int m_Radius;

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

const std::string InitializePathExe::FLAG_NUMBER_OF_REPLICAS = 
  InitializePathExe::FLAG_DELIMITER + "nreplicas";

const std::string InitializePathExe::FLAG_RADIUS = 
  InitializePathExe::FLAG_DELIMITER + "radius";

InitializePathExe::InitializePathExe( int inArgs, char ** iaArgs ) : 
  FreeSurferExecutable( inArgs, iaArgs ) {

  m_EigenVectorFileName = NULL;
  m_SeedVolumeFileName = NULL;
  m_OutputDir = NULL;
  m_NumberOfReplicas = 1;
  m_Radius = 1;

  SetName( "dmri_initialize_path", "find starting path for input to poistats" );  

  SetNextRequiredArgument( FLAG_SEEDS, "seeds", 
    "Seed volume", 
    "seeds.nii", "must specify a seed volume" );

  SetNextRequiredArgument( FLAG_OUTPUT_DIRECTORY, "output", 
    "Output directory", 
    "output/", "must specify an output directory" );

  SetNextOptionalArgument( FLAG_SEED_VALUES, "seednumvalue", 
    "Use <seednumvalue> to define seed region.  If nothing is specified, then all seeds are used.  More than one must be available.  Eg, --seednums 1,2" );

  SetNextOptionalArgument( FLAG_EIGENVECTOR, "eigvec", "Principle eigenvector input" );

  SetNextOptionalArgument( FLAG_NUMBER_OF_REPLICAS, "nreplicas", 
    "Number of times to repeat initialization" );

  SetNextOptionalArgument( FLAG_RADIUS, "radius", "maximum distance the field line initialization can be away" );

  std::string output = "";
  output = output + 
    "InitialPath - volume of initial points" + "\n   " +
    "InitialPath.txt - coordinates of initial path";

  SetOutput( output );
  
  SetVersion( "0.1" );
  SetBugEmail( "martinos-tech@yahoogroups.com" );
    
}

InitializePathExe::~InitializePathExe() {  
}

bool
InitializePathExe::FillArguments() {
  bool isFilled = false;

  try {
    std::string *requiredArguments = GetRequiredArguments();  

    m_SeedVolumeFileName = requiredArguments[0].c_str();
    m_OutputDir = requiredArguments[1].c_str();
    
    isFilled = true;    
  } catch(...) {
    isFilled = false;
  }

  if( isFilled ) {  
    m_SeedValues = m_Parser->GetArgumentIntVector( FLAG_SEED_VALUES.c_str() );
    m_EigenVectorFileName = m_Parser->GetArgument( FLAG_EIGENVECTOR.c_str() );

    // make sure that the number of replicas is 1 or above
    m_NumberOfReplicas = m_Parser->GetArgumentInt( FLAG_NUMBER_OF_REPLICAS.c_str() );
    if( m_NumberOfReplicas < 1 ) {
      m_NumberOfReplicas = 1;
    }
    
    // make sure that the number of replicas is 1 or above
    m_Radius = m_Parser->GetArgumentInt( FLAG_RADIUS.c_str() );
    if( m_Radius < 1 ) {
      m_Radius = 1;
    }
    
  }
  
  return isFilled;
}

void
InitializePathExe::Run() {

  MRI* eigenVectors = NULL;
  
  // create the model for the initializations
  PoistatsModel model;
  model.SetFieldLineRadius( m_Radius );
  model.SetNumberOfControlPoints( 2 );

  // read seed volume
  MRI* seedVolume = FreeSurferExecutable::ReadMRI( m_SeedVolumeFileName );
  model.SetSeedVolume( seedVolume );
  
  // set seed values to connect
  model.SetSeedValues( &m_SeedValues );
  
  // if the eigenvectors were provided, then we'll run the eigenvector intialization
  if( m_EigenVectorFileName != NULL ) {

    // read eigen vectors
    eigenVectors = FreeSurferExecutable::ReadMRI( m_EigenVectorFileName );
    model.SetEigenVectors( eigenVectors );
    
    model.SetInitializePathMode( PoistatsModel::INITIALIZE_EIGEN_VECTOR );
    
  } else {
    model.SetInitializePathMode( PoistatsModel::INITIALIZE_FIELD_LINE );
  }

  InitializePath *initializePath = model.GetPathInitializer();
  
  bool isFirst = true;
  
  for( int i=0; i<m_NumberOfReplicas; i++  ) {
    
    std::cerr << "replica: " << i << std::endl;
    
    // now that all the inputs are set, calculate the initial path
    initializePath->CalculateInitialPath();
    
    // get the resulting pathway
    itk::Array2D< double > *path = initializePath->GetInitialPath();
    
    // save the pathway to the output directory
    const std::string pathFileName( (std::string)m_OutputDir + 
      (std::string)"/InitialPath.txt" );

    if( isFirst ) {
      
      // overwrite previous enteries
      WriteData( pathFileName, path->data_array(), path->rows(), path->cols() );
      isFirst = false;
      
    } else {
      
      // write out the delimiter
      const float delimiter = -999;
      float delimiterPoint[] = { delimiter, delimiter, delimiter };
      WriteDataAppend( pathFileName, delimiterPoint, 1, 3 );
      
      // write the pathway appended
      WriteDataAppend( pathFileName, path->data_array(), path->rows(), path->cols() );
      
    }
        
  }

  if( eigenVectors != NULL ) {
    MRIfree( &eigenVectors );
    eigenVectors = NULL;
  }
  
  if( seedVolume != NULL ) {
    MRIfree( &seedVolume );
    seedVolume = NULL;
  }
    
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
