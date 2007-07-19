#include <iostream>
#include <string>

#include "ui/CommandParser.h"
#include "ui/FreeSurferExecutable.h"

/** This is needed by the freesurfer utils library */
char *Progname;

/**
 * Mask and binarize a path based on the eigen vector.
 */
class MaskExe : public FreeSurferExecutable {

  public:
    
    // string indicating a flag
    static const std::string FLAG_DELIMITER;
    
    // output directory
    static const std::string FLAG_OUTPUT_DIRECTORY;

    // principle eigenvector input
    static const std::string FLAG_EIGENVECTOR;

    // path input
    static const std::string FLAG_PATH;
    
    MaskExe( int inArgs, char ** iaArgs );
    ~MaskExe();

    /**
     * Fills in the arguments and returns true of the arguments can be filled.
     */
    bool FillArguments();    
    
    void Run();
      
  private:

    const char *m_OutputDir;
    const char *m_EigenVectorFileName;
    const char *m_PathFileName;
    
    /**
     * Calculate the vector leading from one point in the path to the next
     */
    void CalculatePathVectors( std::vector< double* > *iPathPoints, 
      std::vector< double* > *oPathVectors );

    /**
     * Mask the path based on the location of the point and the surrounding
     * eigenvectors
     */
    MRI* MaskPath( MRI* iEigenVectors, std::vector< double* > *iPathPoints, 
      std::vector< double* > *iPathVectors );

    /**
     * Accessor for getting the eigenvector at an index
     */
    void GetEigenVector( MRI* iEigenVolume, int *iIndex, double *oVector );
    
    /**
     * Returns true if the eigenvector is a close fit to the path vector.
     */
    bool IsInPath( double *iEigenVector, double *iPathVector );
    

};

const std::string MaskExe::FLAG_DELIMITER = "--";

const std::string MaskExe::FLAG_OUTPUT_DIRECTORY = MaskExe::FLAG_DELIMITER + "o";

const std::string MaskExe::FLAG_EIGENVECTOR = MaskExe::FLAG_DELIMITER + "ev";

const std::string MaskExe::FLAG_PATH = MaskExe::FLAG_DELIMITER + "path";

MaskExe::MaskExe( int inArgs, char ** iaArgs ) : 
  FreeSurferExecutable( inArgs, iaArgs ) {

  m_EigenVectorFileName = NULL;
  m_OutputDir = NULL;

  SetName( "dmri_mask", "masked based on eigenvector and path" );  

  SetNextRequiredArgument( FLAG_OUTPUT_DIRECTORY, "output", 
    "Output directory", 
    "output/", "must specify an output directory" );

  SetNextRequiredArgument( FLAG_EIGENVECTOR, "eigvec", 
    "Principle eigenvector input", 
    "eigvec", "must specify an eigenvector volume" );

  SetNextRequiredArgument( FLAG_PATH, "path", 
    "Path points input", 
    "path", "must specify a path input" );

  std::string output = "";
  output = output + 
    "MaskedPath.mgz - binarized path";

  SetOutput( output );
  
  SetVersion( "0.1" );
  SetBugEmail( "martinos-tech@yahoogroups.com" );
    
}

MaskExe::~MaskExe() {  
}

bool
MaskExe::FillArguments() {
  bool isFilled = false;

  try {
    std::string *requiredArguments = GetRequiredArguments();  

    m_OutputDir = requiredArguments[0].c_str();
    m_EigenVectorFileName = requiredArguments[1].c_str();
    m_PathFileName = requiredArguments[2].c_str();
    
    isFilled = true;    
  } catch(...) {
    isFilled = false;
  }

  if( isFilled ) {      
  }
  
  return isFilled;
}

void
MaskExe::Run() {
  
  // read in eigenvectors
  MRI* eigenVectors = FreeSurferExecutable::ReadMRI( m_EigenVectorFileName );  
  
  // read in points of the path
  std::vector< double* > pathPoints;
  FreeSurferExecutable::ReadData( m_PathFileName, 3, &pathPoints );
    
  // create vectors at the points
  std::vector< double* > pathVectors;
  this->CalculatePathVectors( &pathPoints, &pathVectors );
  
  // mask based on vectors and eigenvectors
  MRI* maskedPath = this->MaskPath( eigenVectors, &pathPoints, &pathVectors );
  
  // save out mask
  const std::string imageFileExtension = ".mgz";
  std::string maskFileName = (std::string)m_OutputDir + 
    (std::string)"/MaskedPath" + imageFileExtension;
  MRIwrite( maskedPath, (char*)maskFileName.c_str() );
  
  // free the vectors
  if( eigenVectors != NULL ) {
    MRIfree( &eigenVectors );
    eigenVectors = NULL;
  }
  
  // free the path points
  for( unsigned int nPoint=0; nPoint<pathPoints.size(); nPoint++ ) {
    double *point = pathPoints[ nPoint ];
    delete []point;
  }
  
  // free the path vector
  for( unsigned int nPoint=0; nPoint<pathVectors.size(); nPoint++ ) {
    double *point = pathVectors[ nPoint ];
    delete []point;
  }
  
  // free the masked path volume
  if( maskedPath != NULL ) {
    MRIfree( &maskedPath );
    maskedPath = NULL;
  }
  
  
}

void 
MaskExe::CalculatePathVectors( std::vector< double* > *iPathPoints, 
  std::vector< double* > *oPathVectors ) {

  const int cCols = 3;
  
  // iterate through all points, but the last
  for( unsigned int nPoint=0; nPoint<iPathPoints->size()-1; nPoint++ ) {
    
    // get the point
    double *point1 = ( *iPathPoints )[ nPoint ];
    
    // get the point ahead of it
    double *point2 = ( *iPathPoints )[ nPoint + 1 ];
    
    // calcualte the difference between the points
    double *difference = new double[ cCols ];
    for( int nCol=0; nCol<cCols; nCol++ ) {
      difference[ nCol ] = point2[ nCol ] - point1[ nCol ];
    }
    
    // get the squared sum
    double squaredSum = 0;
    for( int nCol=0; nCol<cCols; nCol++ ) {
      squaredSum += difference[ nCol ] * difference[ nCol ];
    }
    
    // normalize by the square-rooted squared sum
    double sqrtSum = sqrt( squaredSum );
    for( int nCol=0; nCol<cCols; nCol++ ) {
      difference[ nCol ] /= sqrtSum;
    }
    
    // the difference is the vector that we want to save
    oPathVectors->push_back( difference );
    
  }
    
}

MRI* 
MaskExe::MaskPath( MRI* iEigenVectors, std::vector< double* > *iPathPoints, 
  std::vector< double* > *iPathVectors ) {

  // allocate the mask
  MRI* mask = MRIalloc( iEigenVectors->width, iEigenVectors->height, 
    iEigenVectors->depth, iEigenVectors->type );
  MRIcopyHeader( iEigenVectors, mask );
  
  // set all the pixels in the mask to 0 initially
  for(int z =0 ; z < mask->depth ; z++)
    for(int y = 0 ; y < mask->height ; y++)
      for(int x = 0 ; x < mask->width ; x++)
        MRIFvox( mask, x, y, z) = 0.0f;

  // search radius in voxels around each point that that the current vector will
  // be compared
  // TODO: should I add the radius to the command line?
  const int radius = 1;
  
  // go through all the vectors (which is one less than the path points
  for( unsigned int n=0; n<iPathVectors->size(); n++ ) {
    
    // for a vector associated with a point, search around it and compare the
    // eigenvectors
    double *pathPoint = ( *iPathPoints )[ n ];
    double *pathVector = ( *iPathVectors )[ n ];
    
    // round the point to an int, so it can be an index into the eigen vector 
    // volume
    const int roundedPoint[3] = {
      static_cast< int >( round( pathPoint[ 0 ] ) ),
      static_cast< int >( round( pathPoint[ 1 ] ) ),
      static_cast< int >( round( pathPoint[ 2 ] ) )
    };
    
    // TODO: make sure that the bounds are in the volume
    // calculate the bounds to look in
    const int xBound[2] = { roundedPoint[0] - radius, roundedPoint[0] + radius };
    const int yBound[2] = { roundedPoint[1] - radius, roundedPoint[1] + radius };
    const int zBound[2] = { roundedPoint[2] - radius, roundedPoint[2] + radius };
    
    // look inclusively in the bounds and compare the eigenvectors
    for( int x=xBound[0]; x<=xBound[1]; x++ ) {
      for( int y=yBound[0]; y<=yBound[1]; y++ ) {
        for( int z=zBound[0]; z<=zBound[1]; z++ ) {
          
          // get the eigen vector at the coordinate
          double eigenVector[ 3 ];
          int index[ 3 ] = { x, y, z };
          this->GetEigenVector( iEigenVectors, index, eigenVector );
          
          // compare the path vector to the eigenvector
          if( this->IsInPath( eigenVector, pathVector ) ) {
            MRIFvox( mask, x, y, z) = 1.0;
          }
          
        }
      }
    }
    
  }
  
  return mask;
    
}

void 
MaskExe::GetEigenVector( MRI* iEigenVolume, int *iIndex, double *oVector ) {

  for( int nDim=0; nDim<3; nDim++ ) {
    // get the eigenvector at the current point
    oVector[ nDim ] = MRIgetVoxVal( iEigenVolume, iIndex[ 0 ], iIndex[ 1 ], 
      iIndex[ 2 ], nDim );      
  }
  
  // normalize the vector
  double squaredSum = 0.0;
  for( int nDim=0; nDim<3; nDim++ ) {
    squaredSum += oVector[ nDim ] * oVector[ nDim ];
  }
  
  double sum = sqrt( squaredSum );
  for( int nDim=0; nDim<3; nDim++ ) {
    oVector[ nDim ] /= sum;
  }
  
}

bool 
MaskExe::IsInPath( double *iEigenVector, double *iPathVector ) {
  
  bool isInPath = false;
  
  // take the dot product
  double dotProduct = 0.0;
  for( int nDim=0; nDim<3; nDim++ ) {
    dotProduct += iEigenVector[ nDim ] * iPathVector[ nDim ];
  }
    
  // this means it's greater than 90 degrees
  if( dotProduct < 0 ) {
    
    // recalculate with flipped the eigenvector
    dotProduct = 0.0;
    for( int nDim=0; nDim<3; nDim++ ) {
      dotProduct += -iEigenVector[ nDim ] * iPathVector[ nDim ];
    }
    
  }
      
  const double angleBetweenVectors = acos( dotProduct );

  const double threshold = PI / 4.0;  
  
  // if the dot product is too big, we're not in the path
  // TODO: figure out what a good number is...
  if( angleBetweenVectors < threshold ) {
    isInPath = true;
  }  
  
  return isInPath;
}

int main( int argc, char ** argv ) {
  
  // this is needed by the freesurfer utils library
  Progname = argv[0];
  
  FreeSurferExecutable *exe = new MaskExe( argc, argv );
  
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
