#ifndef _itkPoistatsFilter_txx
#define _itkPoistatsFilter_txx

#include <itkBSplineInterpolateImageFunction.h>
#include <itkPointSet.h>

// for calculating the elapsed time
#include <time.h>

#include <vnl/vnl_diag_matrix.h>
#include <vnl/vnl_cross.h>

#include "itkPoistatsFilter.h"

namespace itk
{

// I put two opening and closing braces around each tuple in order to silence 
// a compile warning
template <class TInputImage, class TOutputImage>
const double
PoistatsFilter<TInputImage, TOutputImage>
::NO_ZERO_SHELL_252[][3] = {
    {-0.99419,0,-0.107642},
    {-0.985146,-0.130799,0.111264},
    {-0.985146,0.130799,0.111264},
    {-0.963803,-0.245346,-0.104352},
    {-0.963803,0.245346,-0.104352},
    {-0.948206,-0.115895,-0.295759},
    {-0.948206,0.115895,-0.295759},
    {-0.942646,0,0.333793},
    {-0.923944,-0.36802,0.104352},
    {-0.923944,0.36802,0.104352},
    {-0.911958,-0.253082,0.322926},
    {-0.911958,0.253082,0.322926},
    {-0.894427,0,-0.447214},
    {-0.886548,-0.331288,-0.322926},
    {-0.886548,0.331288,-0.322926},
    {-0.873881,-0.473235,-0.111264},
    {-0.873881,0.473235,-0.111264},
    {-0.84962,-0.187522,-0.492931},
    {-0.84962,0.187522,-0.492931},
    {-0.844226,-0.122673,0.52176},
    {-0.844226,0.122673,0.52176},
    {-0.835236,-0.463581,0.295759},
    {-0.835236,0.463581,0.295759},
    {-0.804316,-0.58437,0.107642},
    {-0.804316,0.58437,0.107642},
    {-0.79758,-0.347685,0.492931},
    {-0.79758,0.347685,0.492931},
    {-0.78869,0,-0.614791},
    {-0.762617,-0.554074,-0.333793},
    {-0.762617,0.554074,-0.333793},
    {-0.755099,-0.396979,-0.52176},
    {-0.755099,0.396979,-0.52176},
    {-0.723607,-0.525731,0.447214},
    {-0.723607,0.525731,0.447214},
    {-0.720118,-0.684873,-0.111264},
    {-0.720118,0,0.693852},
    {-0.720118,0.684873,-0.111264},
    {-0.712379,-0.204747,-0.671263},
    {-0.712379,0.204747,-0.671263},
    {-0.698993,-0.651103,0.295759},
    {-0.698993,0.651103,0.295759},
    {-0.696674,-0.253082,0.671263},
    {-0.696674,0.253082,0.671263},
    {-0.638064,-0.463581,0.614791},
    {-0.638064,0.463581,0.614791},
    {-0.635522,-0.764998,0.104352},
    {-0.635522,0.764998,0.104352},
    {-0.626112,0,-0.779733},
    {-0.610888,-0.595468,-0.52176},
    {-0.610888,0.595468,-0.52176},
    {-0.589032,-0.740783,-0.322926},
    {-0.589032,0.740783,-0.322926},
    {-0.582588,-0.423275,-0.693852},
    {-0.582588,0.423275,-0.693852},
    {-0.577134,-0.651103,0.492931},
    {-0.577134,0.651103,0.492931},
    {-0.540088,-0.130799,0.831382},
    {-0.540088,0.130799,0.831382},
    {-0.53117,-0.840815,-0.104352},
    {-0.53117,0.840815,-0.104352},
    {-0.522506,-0.789117,0.322926},
    {-0.522506,0.789117,0.322926},
    {-0.513822,-0.211637,-0.831382},
    {-0.513822,0.211637,-0.831382},
    {-0.506536,-0.36802,0.779733},
    {-0.506536,0.36802,0.779733},
    {-0.455979,-0.58437,0.671263},
    {-0.455979,0.58437,0.671263},
    {-0.440891,-0.750089,-0.492931},
    {-0.440891,0.750089,-0.492931},
    {-0.428824,-0.89651,0.111264},
    {-0.428824,0.89651,0.111264},
    {-0.417408,0,-0.908719},
    {-0.414864,-0.614242,-0.671263},
    {-0.414864,0.614242,-0.671263},
    {-0.403235,-0.865984,-0.295759},
    {-0.403235,0.865984,-0.295759},
    {-0.377549,-0.764998,0.52176},
    {-0.377549,0.764998,0.52176},
    {-0.360059,-0.423275,-0.831382},
    {-0.360059,0.423275,-0.831382},
    {-0.348337,0,0.937369},
    {-0.33769,-0.245346,0.908719},
    {-0.33769,0.245346,0.908719},
    {-0.307222,-0.945531,-0.107642},
    {-0.307222,0.945531,-0.107642},
    {-0.291294,-0.89651,0.333793},
    {-0.291294,-0.473235,0.831382},
    {-0.291294,0.473235,0.831382},
    {-0.291294,0.89651,0.333793},
    {-0.281811,-0.204747,-0.937369},
    {-0.281811,0.204747,-0.937369},
    {-0.276393,-0.850651,-0.447214},
    {-0.276393,0.850651,-0.447214},
    {-0.243719,-0.750089,-0.614791},
    {-0.243719,0.750089,-0.614791},
    {-0.222529,-0.684873,0.693852},
    {-0.222529,0.684873,0.693852},
    {-0.197173,0,-0.980369},
    {-0.193479,-0.595468,-0.779733},
    {-0.193479,0.595468,-0.779733},
    {-0.182789,-0.937611,-0.295759},
    {-0.182789,0.937611,-0.295759},
    {-0.180029,-0.977348,0.111264},
    {-0.180029,0.977348,0.111264},
    {-0.159516,-0.115895,0.980369},
    {-0.159516,0.115895,0.980369},
    {-0.144211,-0.840815,0.52176},
    {-0.144211,0.840815,0.52176},
    {-0.128986,-0.396979,-0.908719},
    {-0.128986,0.396979,-0.908719},
    {-0.107642,-0.331288,0.937369},
    {-0.107642,0.331288,0.937369},
    {-0.0842027,-0.865984,-0.492931},
    {-0.0842027,0.865984,-0.492931},
    {-0.0644931,-0.992447,-0.104352},
    {-0.0644931,0.992447,-0.104352},
    {-0.0609297,-0.187522,-0.980369},
    {-0.0609297,0.187522,-0.980369},
    {-0.0424992,-0.554074,0.831382},
    {-0.0424992,0.554074,0.831382},
    {-0.0411156,-0.945531,0.322926},
    {-0.0411156,0.945531,0.322926},
    {-0.0254109,-0.740783,-0.671263},
    {-0.0254109,0.740783,-0.671263},
    {0,0,-1.},
    {0,0,1.},
    {0.0254109,-0.740783,0.671263},
    {0.0254109,0.740783,0.671263},
    {0.0411156,-0.945531,-0.322926},
    {0.0411156,0.945531,-0.322926},
    {0.0424992,-0.554074,-0.831382},
    {0.0424992,0.554074,-0.831382},
    {0.0609297,-0.187522,0.980369},
    {0.0609297,0.187522,0.980369},
    {0.0644931,-0.992447,0.104352},
    {0.0644931,0.992447,0.104352},
    {0.0842027,-0.865984,0.492931},
    {0.0842027,0.865984,0.492931},
    {0.107642,-0.331288,-0.937369},
    {0.107642,0.331288,-0.937369},
    {0.128986,-0.396979,0.908719},
    {0.128986,0.396979,0.908719},
    {0.144211,-0.840815,-0.52176},
    {0.144211,0.840815,-0.52176},
    {0.159516,-0.115895,-0.980369},
    {0.159516,0.115895,-0.980369},
    {0.180029,-0.977348,-0.111264},
    {0.180029,0.977348,-0.111264},
    {0.182789,-0.937611,0.295759},
    {0.182789,0.937611,0.295759},
    {0.193479,-0.595468,0.779733},
    {0.193479,0.595468,0.779733},
    {0.197173,0,0.980369},
    {0.222529,-0.684873,-0.693852},
    {0.222529,0.684873,-0.693852},
    {0.243719,-0.750089,0.614791},
    {0.243719,0.750089,0.614791},
    {0.276393,-0.850651,0.447214},
    {0.276393,0.850651,0.447214},
    {0.281811,-0.204747,0.937369},
    {0.281811,0.204747,0.937369},
    {0.291294,-0.89651,-0.333793},
    {0.291294,-0.473235,-0.831382},
    {0.291294,0.473235,-0.831382},
    {0.291294,0.89651,-0.333793},
    {0.307222,-0.945531,0.107642},
    {0.307222,0.945531,0.107642},
    {0.33769,-0.245346,-0.908719},
    {0.33769,0.245346,-0.908719},
    {0.348337,0,-0.937369},
    {0.360059,-0.423275,0.831382},
    {0.360059,0.423275,0.831382},
    {0.377549,-0.764998,-0.52176},
    {0.377549,0.764998,-0.52176},
    {0.403235,-0.865984,0.295759},
    {0.403235,0.865984,0.295759},
    {0.414864,-0.614242,0.671263},
    {0.414864,0.614242,0.671263},
    {0.417408,0,0.908719},
    {0.428824,-0.89651,-0.111264},
    {0.428824,0.89651,-0.111264},
    {0.440891,-0.750089,0.492931},
    {0.440891,0.750089,0.492931},
    {0.455979,-0.58437,-0.671263},
    {0.455979,0.58437,-0.671263},
    {0.506536,-0.36802,-0.779733},
    {0.506536,0.36802,-0.779733},
    {0.513822,-0.211637,0.831382},
    {0.513822,0.211637,0.831382},
    {0.522506,-0.789117,-0.322926},
    {0.522506,0.789117,-0.322926},
    {0.53117,-0.840815,0.104352},
    {0.53117,0.840815,0.104352},
    {0.540088,-0.130799,-0.831382},
    {0.540088,0.130799,-0.831382},
    {0.577134,-0.651103,-0.492931},
    {0.577134,0.651103,-0.492931},
    {0.582588,-0.423275,0.693852},
    {0.582588,0.423275,0.693852},
    {0.589032,-0.740783,0.322926},
    {0.589032,0.740783,0.322926},
    {0.610888,-0.595468,0.52176},
    {0.610888,0.595468,0.52176},
    {0.626112,0,0.779733},
    {0.635522,-0.764998,-0.104352},
    {0.635522,0.764998,-0.104352},
    {0.638064,-0.463581,-0.614791},
    {0.638064,0.463581,-0.614791},
    {0.696674,-0.253082,-0.671263},
    {0.696674,0.253082,-0.671263},
    {0.698993,-0.651103,-0.295759},
    {0.698993,0.651103,-0.295759},
    {0.712379,-0.204747,0.671263},
    {0.712379,0.204747,0.671263},
    {0.720118,-0.684873,0.111264},
    {0.720118,0,-0.693852},
    {0.720118,0.684873,0.111264},
    {0.723607,-0.525731,-0.447214},
    {0.723607,0.525731,-0.447214},
    {0.755099,-0.396979,0.52176},
    {0.755099,0.396979,0.52176},
    {0.762617,-0.554074,0.333793},
    {0.762617,0.554074,0.333793},
    {0.78869,0,0.614791},
    {0.79758,-0.347685,-0.492931},
    {0.79758,0.347685,-0.492931},
    {0.804316,-0.58437,-0.107642},
    {0.804316,0.58437,-0.107642},
    {0.835236,-0.463581,-0.295759},
    {0.835236,0.463581,-0.295759},
    {0.844226,-0.122673,-0.52176},
    {0.844226,0.122673,-0.52176},
    {0.84962,-0.187522,0.492931},
    {0.84962,0.187522,0.492931},
    {0.873881,-0.473235,0.111264},
    {0.873881,0.473235,0.111264},
    {0.886548,-0.331288,0.322926},
    {0.886548,0.331288,0.322926},
    {0.894427,0,0.447214},
    {0.911958,-0.253082,-0.322926},
    {0.911958,0.253082,-0.322926},
    {0.923944,-0.36802,-0.104352},
    {0.923944,0.36802,-0.104352},
    {0.942646,0,-0.333793},
    {0.948206,-0.115895,0.295759},
    {0.948206,0.115895,0.295759},
    {0.963803,-0.245346,0.104352},
    {0.963803,0.245346,0.104352},
    {0.985146,-0.130799,-0.111264},
    {0.985146,0.130799,-0.111264},
    {0.99419,0,0.107642}
  };
  

/**
 * Constructor
 */
template <class TInputImage, class TOutputImage>
PoistatsFilter<TInputImage, TOutputImage>
::PoistatsFilter()
{

  std::srand( ( unsigned ) time( 0 ) );
  const long seed = std::rand();

  this->m_PoistatsModel = new PoistatsModel( seed );  

  this->SetNumberOfRequiredOutputs( 2 );
      
  this->SetNumberOfSamplePoints( DEFAULT_NUMBER_OF_SAMPLE_POINTS );
  
  this->SetMaxLull( DEFAULT_MAX_LULL );
  this->SetMaxTime( DEFAULT_MAX_TIME );
  
  this->SetCoolFactor( DEFAULT_COOL_FACTOR );
  
  this->SetNumberOfDirections( DEFAULT_NUMBER_OF_DIRECTIONS );

  this->m_OdfLookUpTable = NULL;

  this->m_Polarity = MatrixType( 3, 3 );
  this->m_Polarity.fill( 1 );
  
  m_Replicas = new PoistatsReplicas( this->m_PoistatsModel, 
    DEFAULT_NUMBER_OF_REPLICAS );
    
  m_Replicas->SetNumberOfSteps( DEFAULT_NUMBER_OF_STEPS );
  
  const double smallValue = 1e-200;
  m_InvalidOdf = ArrayType( GetNumberOfDirections() );
  m_InvalidOdf.fill( smallValue );
  
  this->SetReplicaExchangeProbability( DEFAULT_REPLICA_EXCHANGE_PROBABILITY );

  this->SetSigmaTimeConstant( DEFAULT_SIGMA_TIME_CONSTANT );

  this->SetPointsToImageGamma( DEFAULT_POINTS_TO_IMAGE_GAMMA );
  
}

template <class TInputImage, class TOutputImage>
PoistatsFilter<TInputImage, TOutputImage>
::~PoistatsFilter() {
  delete this->m_Replicas;
  delete m_PoistatsModel;  
}

template <class TInputImage, class TOutputImage>
int
PoistatsFilter<TInputImage, TOutputImage>
::GetNumberOfReplicas() {
  return this->m_Replicas->GetNumberOfReplicas();
}

template <class TInputImage, class TOutputImage>
void
PoistatsFilter<TInputImage, TOutputImage>
::SetNumberOfReplicas( const int nReplicas ) {
  this->m_Replicas->SetNumberOfReplicas( nReplicas );
}


template <class TInputImage, class TOutputImage>
typename PoistatsFilter<TInputImage, TOutputImage>::OdfLookUpTablePointer
PoistatsFilter<TInputImage, TOutputImage>
::GetOdfLookUpTable() {

  if( !this->m_OdfLookUpTable ) {
  
    this->m_OdfLookUpTable = OdfLookUpTableType::New();
          
    RegionType dtiRegion = this->GetInput()->GetLargestPossibleRegion();
    OdfLookUpRegionType odfRegion;
    
    double odfOrigin[ OdfLookUpRegionType::GetImageDimension() ];
    OdfLookUpIndexType odfStart;

    OdfLookUpSizeType odfSize;
    for( unsigned int cDim=0; cDim<OdfLookUpRegionType::GetImageDimension(); cDim++ ) {    
      odfSize[ cDim ] = dtiRegion.GetSize()[ cDim ];
      odfOrigin[ cDim ] = this->GetInput()->GetOrigin()[ cDim ];
      odfStart[ cDim ] = 0;
    }
        
    odfRegion.SetSize( odfSize );    
    odfRegion.SetIndex( odfStart );
                    
    this->m_OdfLookUpTable->SetRegions( odfRegion );
    this->m_OdfLookUpTable->SetOrigin( odfOrigin );
    
    this->m_OdfLookUpTable->Allocate();

    // initially set all the indices to be invalid
    this->m_OdfLookUpTable->FillBuffer( INVALID_INDEX );
  }

  return this->m_OdfLookUpTable;
}

template <class TInputImage, class TOutputImage>
void
PoistatsFilter<TInputImage, TOutputImage>
::AllocateOutputImage( const int outputIndex, OutputImagePointer image ) {
  
  InputImageConstPointer inputImage = this->GetInput();
  
  RegionType dtiRegion = inputImage->GetLargestPossibleRegion();
  OutputSizeType outputSize;
  double outputOrigin[ OutputRegionType::GetImageDimension() ];
  OutputIndexType outputStart;
  OutputSpacingType outputSpacing;

  for( unsigned int cDim=0; cDim<OutputRegionType::GetImageDimension(); cDim++ ) {    
    outputSize[ cDim ] = dtiRegion.GetSize()[ cDim ];
    outputOrigin[ cDim ] = this->GetInput()->GetOrigin()[ cDim ];
    outputStart[ cDim ] = 0;
    outputSpacing[ cDim ] = this->GetInput()->GetSpacing()[ cDim ];
  }
      
  OutputRegionType outputRegion;
  outputRegion.SetSize( outputSize );    
  outputRegion.SetIndex( outputStart );
  
  image->SetRegions( outputRegion );
  image->SetOrigin( outputOrigin );
  image->SetSpacing( outputSpacing );
  
  // set the cosine directions
  image->SetDirection( this->GetInput()->GetDirection() );
  
  image->Allocate();

  image->FillBuffer( 0.0 );

  this->SetNthOutput( outputIndex, image );
  
}

template <class TInputImage, class TOutputImage>
void
PoistatsFilter< TInputImage, TOutputImage >
::GetPositiveMinimumInt( MatrixPointer points, const int radius, 
  VoxelIndexPointer minimum ) {

  const int rowFloor = 0;
  
  for( unsigned int cColumn=0; cColumn<points->cols(); cColumn++ ) {

    vnl_vector< double > column = points->get_column( cColumn );
    int minimumOfRow = static_cast< int >( floor( column.min_value() ) );
    minimumOfRow -= radius;
    
    if( minimumOfRow < rowFloor ) {
      minimumOfRow = rowFloor;
    }
    
    ( *minimum )[ cColumn ] = minimumOfRow;
        
  }
  
}

template <class TInputImage, class TOutputImage>
void
PoistatsFilter< TInputImage, TOutputImage >
::GetPositiveMaximumInt( MatrixPointer points, const int radius, 
  VoxelIndexPointer maximum, OutputRegionType pointsRegion ) {

  for( unsigned int cColumn=0; cColumn<points->cols(); cColumn++ ) {
    
    vnl_vector< double > column = points->get_column( cColumn );
    int maximumOfRow = static_cast< int >( ceil( column.max_value() ) );
    maximumOfRow += radius;
    
    const int rowCeiling = pointsRegion.GetSize()[ cColumn ] - 1;
    if( maximumOfRow > rowCeiling ) {
      maximumOfRow = rowCeiling;
    }
    
    ( *maximum )[ cColumn ] = maximumOfRow;
    
  }

}



template <class TInputImage, class TOutputImage>
void
PoistatsFilter< TInputImage, TOutputImage >
::ConvertPointsToImage( MatrixPointer points, OutputImagePointer image ) {

  /* MATLAB:  
  mn = floor(min(pnts)-3);
  mn = max(1,mn);
  
  mx = ceil(max(pnts)+3);
  mx = min(volsize, mx);
  
  [X Y Z] = meshgrid(mn(1):mx(1), mn(2):mx(2), mn(3):mx(3));
  X = X(:); Y = Y(:); Z = Z(:);
  box = [X Y Z];
  wts = sum(exp(-gamma*distmat(box, pnts).^2),2);
  wts = wts/sum(wts(:));
  idx = int64(sub2ind(volsize, X, Y, Z));
  
  vol = zeros(volsize);
  vol(idx) = vol(idx) + wts;
  */

  const int radius = 3;

  VoxelIndexType minimumPoint( 3 );
  GetPositiveMinimumInt( points, radius, &minimumPoint );
  
  VoxelIndexType maximumPoint( 3 );
  GetPositiveMaximumInt( points, radius, &maximumPoint, 
    image->GetLargestPossibleRegion() );

  const int y = 0;
  const int x = 1;
  const int z = 2;
    
  const int columnStart = minimumPoint[ y ];
  const int columnEnd = maximumPoint[ y ];

  const int rowStart = minimumPoint[ x ];
  const int rowEnd = maximumPoint[ x ];
  
  const int sliceStart = minimumPoint[ z ];
  const int sliceEnd = maximumPoint[ z ];

  const double gamma = this->GetPointsToImageGamma();

  itk::Array< double > box( points->cols() );
  OutputIndexType pixelIndex;

  double totalWeights = 0.0;

  // these for loops are inclusive of the ending indices

  for( int cCol=columnStart; cCol<=columnEnd; cCol++ ) {
    for( int cRow=rowStart; cRow<=rowEnd; cRow++ ) {      
      for( int cSlice=sliceStart; cSlice<=sliceEnd; cSlice++ ) {

        box[ y ] = static_cast< double >( cCol );
        box[ x ] = static_cast< double >( cRow );
        box[ z ] = static_cast< double >( cSlice );

        ArrayType distanceMatrix( points->rows() );

        for( unsigned int cPoint=0; cPoint<points->rows(); cPoint++ ) {

          vnl_vector< double > point = points->get_row( cPoint );

          const double distance = GetDistance( &point, &box );
          distanceMatrix[ cPoint ] = exp( -gamma * distance * distance );

        }
        
        double sum = distanceMatrix.sum();

        pixelIndex[ 0 ] = cCol; 
        pixelIndex[ 1 ] = cRow; 
        pixelIndex[ 2 ] = cSlice;
        
        image->SetPixel( pixelIndex, sum );
        totalWeights += sum;      
      }
    }  
  }

  // divide each element by the totalWeights
  typedef itk::ImageRegionIterator< OutputImageType > IteratorType;
  IteratorType it( image, image->GetLargestPossibleRegion() );
  
  for( it = it.Begin(); !it.IsAtEnd(); ++it ) {
    double pixelValue = it.Value();
    if( pixelValue != 0.0 ) {
      pixelValue /= totalWeights;
      it.Set( pixelValue );
    } 
  }
  
}

/**
 * This is called when Update() is called.  This is were the work of finding
 * the most probable path is done and were all the outputs are generated.
 */
template <class TInputImage, class TOutputImage>
void
PoistatsFilter< TInputImage, TOutputImage >
::GenerateData() {

  this->InvokeEvent( StartEvent() );
    
  itkDebugMacro( << "parsing seed volume" );  
  // extracts the start and end regions of interest
  try {
    this->ParseSeedVolume();
  } catch( itk::ExceptionObject & excp ) {
    throw excp;
  }
  
  // if a mask is being used, make it the union of the seed regions and the mask
  // itself, in case the mask doesn't include the seed regions in it
  this->TakeUnionOfMaskAndSeeds();

  itkDebugMacro( << "initializing paths" );
  // initializes all the paths of the replicas, connects start and end regions
  this->InitPaths();
  
  // creates odfs throughout tensor volume
  this->ConstructOdfList();

  // initialize temperatures to be regularily spaced between 0.05 and 0.1
  const double temperatureFloor = 0.05;
  const double temperatureCeiling = 0.1;
  m_Replicas->SpaceTemperaturesEvenly( temperatureFloor, temperatureCeiling );
  
  this->m_Replicas->FillCurrentMeanEnergies( 0.0 );

  const double maxDouble = std::numeric_limits< double >::max();
  
  this->m_Replicas->FillPreviousMeanEnergies( maxDouble );
  
  this->SetGlobalMinEnergy( maxDouble );

  const int numberOfSpatialDimensions = 3;  
  MatrixType finalBestPath( this->m_Replicas->GetNumberOfSteps(), 
    numberOfSpatialDimensions );
  
  const double lullEnergyDifferenceThreshold = 0.0005;
  
  const bool isMoreThanOneReplica = this->m_Replicas->GetNumberOfReplicas() > 1;

  // this iterates until a minimum is found or we iterate too much      
  for( this->m_CurrentIteration=1, this->m_CurrentLull=0;
       
       this->m_CurrentIteration < this->GetMaxTime() && 
       this->m_CurrentLull < this->GetMaxLull();
       
       this->m_CurrentIteration++ )
    {

    // start the clock
    clock_t startClock = clock();

    // the amount in terms of voxels that the path can perturb
    const double sigma = static_cast< double >( this->GetInitialSigma() ) * 
      exp( -static_cast< double >( m_CurrentIteration ) / 
           this->GetSigmaTimeConstant() );
                   
    // reset the number of exchanges that occured...if you're keeping track
    this->SetExchanges( 0 );

    // now go through all the replicas and wiggle them around    
    for( int cReplica=0; cReplica<this->GetNumberOfReplicas(); cReplica++ ) {
    
      //if time > 1, prevpath{i} = trialpath{i}; end;
      const bool isFirst = m_CurrentIteration == 1;
      if( !isFirst ) {
        this->m_Replicas->CopyCurrentToPreviousTrialPath( cReplica );
      }

      // get the low resolution path for this replica
      MatrixPointer currentBasePath = this->m_Replicas->GetBasePath( cReplica );
      MatrixType lowTrialPath( 
        currentBasePath->rows(), currentBasePath->cols() );
      m_Replicas->GetPerturbedBasePath( cReplica, &lowTrialPath, sigma, 
        this->GetStartSeeds(), this->GetEndSeeds() );
        
      // MATLAB: trialpath{i} = rethreadpath(lowtrialpath, steps);
      this->m_Replicas->PerturbCurrentTrialPath( cReplica, &lowTrialPath, 
        this->GetNumberOfSteps() );
      MatrixPointer perturbedTrialPath = 
        this->m_Replicas->GetCurrentTrialPath( cReplica );

      // MATLAB: rpath = round(trialpath{i});
      itk::Array2D< int > roundedPath( perturbedTrialPath->rows(),  
                                       perturbedTrialPath->cols() );
      
      // we want to obtain an index into our image, so round the path
      this->RoundPath( &roundedPath, perturbedTrialPath );
      
      ArrayPointer odfs[ this->GetNumberOfSteps() ];
      this->GetOdfsAtPoints( odfs, &roundedPath );
      
      /* MATLAB: 
        % calculate path energy      
        energy(i) = odfpathenergy(trialpath{i}, odfs, geo);
      */
      const double meanPathEnergy = 
        this->CalculateOdfPathEnergy( perturbedTrialPath, odfs, NULL );

      this->m_Replicas->SetCurrentMeanEnergy( cReplica, meanPathEnergy );

      /* MATLAB:                   
        % check for Metropolis-Hastings update
        Delta = (energy(i)-energyprev(i))/temp(i);
        updateprobability = min(1, exp(-Delta));
        if rand(1) <= updateprobability;
      */   
      // update the energy, and always set it initially to the first replica
      if( this->m_Replicas->ShouldUpdateEnergy( cReplica ) || 
        ( cReplica == 0 &&  isFirst ) ) {
      
        /* MATLAB:
          basepath{i} = lowtrialpath; 
          bestpath{i} = trialpath{i};
        */
        this->m_Replicas->FoundBestPath( cReplica, &lowTrialPath );

        /* MATLAB:                
          if energy(i) < globalminenergy
            globalbestpath = trialpath{i};          
            globalminenergy = energy(i);
          end
        */
        if( this->m_Replicas->GetCurrentMeanEnergy( cReplica ) < 
          this->GetGlobalMinEnergy() ) {

          PoistatsReplica::CopyPath( perturbedTrialPath, &finalBestPath );

          this->SetGlobalMinEnergy( 
            this->m_Replicas->GetCurrentMeanEnergy( cReplica ) );
          
        }

      } else {
        // MATLAB: energy(i) = energyprev(i);
        this->m_Replicas->ResetCurrentToPreviousEnergy( cReplica );
      }
      
      // MATLAB: if  time > 1 & rand(1) < replicaexchprob
      const bool shouldExchange = !isFirst && isMoreThanOneReplica &&
        ( this->m_PoistatsModel->GetRandomNumber() < 
          this->GetReplicaExchangeProbability() );

      if( shouldExchange ) {
      
        // get a random replica and the next replica with the next highest mean
        int randomReplicaIndex;
        int followingRandomReplicaIndex;
        this->m_Replicas->GetRandomSortedFirstSecondReplicas( 
          randomReplicaIndex, followingRandomReplicaIndex );
          
        const double probabilityExchange = 
          this->m_Replicas->CalculateProbablityExchange( randomReplicaIndex, 
          followingRandomReplicaIndex );

        if( this->m_PoistatsModel->GetRandomNumber() <= probabilityExchange ) {

          this->m_Replicas->ExchangeTemperatures( randomReplicaIndex, 
            followingRandomReplicaIndex );
          this->m_Exchanges++;
            
        }
        
      }
          
    } // end of peturbing every replica for a single iteration
    
    // MATLAB:denergy = (mean(energy)-mean(energyprev))/mean(energy);        
    m_CurrentEnergyDifference = 
      this->m_Replicas->GetNormalizedMeanCurrentPreviousEnergiesDifference();
      
    /* MATLAB:
      if abs(denergy) < lulldenergy
        lull = lull + 1;
      else
        lull = 0;
      end
    */
    if( fabs( m_CurrentEnergyDifference ) < lullEnergyDifferenceThreshold ) {
      m_CurrentLull++;
    } else {
      m_CurrentLull = 0;
    }
        
    // MATLAB: energyprev = energy;
    this->m_Replicas->CopyCurrentToPreviousEnergies();
        
    // MATLAB: temp = coolfactor*temp;     
    this->m_Replicas->CoolTemperatures( this->GetCoolFactor() );

    clock_t endClock = clock();
    const double elapsedTime = 
      static_cast< double >( endClock - startClock ) / CLOCKS_PER_SEC;
    this->SetElapsedTime( elapsedTime );
    
    // let our observers know that we've iterated
    this->InvokeEvent( IterationEvent() );            
    
  } // we've found an optimal path

  // our user has specified the number of sample point to use
  MatrixPointer rethreadedFinalPath = this->m_Replicas->RethreadPath( 
    &finalBestPath, this->GetNumberOfSamplePoints() );
    
  // this is the path that we'll output to the user
  this->SetFinalPath( *rethreadedFinalPath );
  finalBestPath = this->GetFinalPath();
  
  delete rethreadedFinalPath;
  rethreadedFinalPath = NULL;
  
  // this is a volume of our best path
  this->InvokeEvent( GenerateOptimalPathDensitiesEvent() );
  OutputImagePointer optimalPathDensity = OutputImageType::New();
  AllocateOutputImage( OPTIMAL_PATH_DENSITY_OUTPUT, optimalPathDensity );
  this->ConvertPointsToImage( &finalBestPath, optimalPathDensity );  

  // this is an aggregate of all of the best paths found in all the replicas
  this->InvokeEvent( GenerateBestReplicaPathDensitiesStartEvent() );
  OutputImagePointer pathDensities = OutputImageType::New();
  AllocateOutputImage( PATH_DENSITY_OUTPUT, pathDensities );
  MatrixListType bestPaths = m_Replicas->GetBestTrialPaths();
  this->GetAggregateReplicaDensities( bestPaths, pathDensities );
  this->InvokeEvent( GenerateBestReplicaPathDensitiesEndEvent() );
  
  // this is the probabilities associated with the best path at each sample 
  // voxel
  this->InvokeEvent( GenerateFinalPathProbabilitiesEvent() );
  this->CalculateFinalPathProbabilities();
  
//  std::cerr << "calculate best path probabilities for each replica..." << std::endl;
//  this->CalculateBestPathProbabilities();
  
  this->InvokeEvent( EndEvent() );

}

template <class TInputImage, class TOutputImage>
double 
PoistatsFilter<TInputImage, TOutputImage>
::GetCurrentMeanOfEnergies() const {
  return this->m_Replicas->GetCurrentMeanOfEnergies();
}

template <class TInputImage, class TOutputImage>
double 
PoistatsFilter<TInputImage, TOutputImage>
::GetCurrentMinOfEnergies() const {
  return this->m_Replicas->GetMinimumCurrentEnergy();
}

template <class TInputImage, class TOutputImage>
void 
PoistatsFilter<TInputImage, TOutputImage>
::SetMghEigenVectors( MRI *vectors ) {
  this->m_PoistatsModel->SetEigenVectors( vectors );
  this->m_PoistatsModel->SetUsingPathInitialization( true );
}
  
template <class TInputImage, class TOutputImage>
void 
PoistatsFilter<TInputImage, TOutputImage>
::SetMghSeeds( MRI* seeds ) {
  this->m_PoistatsModel->SetSeedVolume( seeds );
}


template <class TInputImage, class TOutputImage>
typename PoistatsFilter<TInputImage, TOutputImage>::MatrixType
PoistatsFilter<TInputImage, TOutputImage>
::GetBestPathsProbabilities() {

// TODO: we're using the number of steps rather than the number of sample points
  MatrixType probabilities( this->m_Replicas->GetNumberOfReplicas(), 
                            this->GetNumberOfSteps() );
//                            this->GetNumberOfSamplePoints() );

  this->m_Replicas->GetBestTrialPathsProbabilities( &probabilities );
  
  return probabilities;

}

/**
 * Calculate the probabilities of the best path that was found.
 */
template <class TInputImage, class TOutputImage>
void
PoistatsFilter<TInputImage, TOutputImage>
::CalculateBestPathProbabilities() {

  for( int cReplica = 0; 
    cReplica<this->m_Replicas->GetNumberOfReplicas(); cReplica++ ) {

    MatrixPointer replicaPath = this->m_Replicas->GetBestTrialPath( cReplica );
    
    ArrayType probabilities( replicaPath->rows() );
    this->GetPathProbabilities( replicaPath, &probabilities );
    
    this->m_Replicas->SetBestTrialPathProbabilities( cReplica, &probabilities );
  }  
}

template <class TInputImage, class TOutputImage>
void
PoistatsFilter<TInputImage, TOutputImage>
::ConstructOdfList() {
  /* MATLAB:

   fprintf('Constructing ODF list ... \n');
   flatd = reshape(dtensortmp, [slices*rows*cols 9]);
   odflist = zeros(length(maskidx),size(geo,1));
   for i = 1:length(maskidx)
     % every 10,000th element, print something
     if ~mod(i, 1e4)
         fprintf('%d ',i); 
     end;
     d = reshape(flatd(maskidx(i),:), [3 3]);
     d = d(order, order) .* polarity;
     d = R*d*invR;
     odflist(i,:) = tensor2odf(d, geo, 'angulargaussian')';
   end
   fprintf(' done (%g)\n',toc);
  */

  InputImageConstPointer inputImage = this->GetInput();
    
  SizeType imageSize = inputImage->GetLargestPossibleRegion().GetSize();
//  const int nColumns = imageSize[ 0 ];
//  const int nRows = imageSize[ 1 ];
//  const int nSlices = imageSize[ 2 ];
  
  int cOdfs = 0;
  
  MatrixType polarity = this->GetPolarity();

  /* MATLAB:   
    % rotation matrix from magnet to slice frame
    % only works for oblique axial
    sliceup = [normal_s normal_a 0];
    R = rot3u2v([1 0 0], sliceup);
    invR = inv(R);  
  */
  MatrixType rotationMatrix( 3, 3 );
  GetMagnetToSliceFrameRotation( &rotationMatrix );
  itkDebugMacro( << "magnet to slice frame rotation:" );
  itkDebugMacro( << rotationMatrix );
  
  const MatrixType inverseRotationMatrix =
    vnl_matrix_inverse< double >( rotationMatrix ).inverse();
  const bool isRotationIdentity = rotationMatrix.is_identity();   
     
  OdfLookUpTablePointer odfLookUpTable = this->GetOdfLookUpTable();
  
  this->InvokeEvent( PoistatsOdfCalculationStartEvent() );

  // start the timer
  clock_t startClock = clock();

  itk::ImageRegionIterator< OdfLookUpTableType > odfLookUpTableIt( 
    odfLookUpTable, odfLookUpTable->GetLargestPossibleRegion() );  

  itk::ImageRegionConstIterator< InputImageType > inputImageIt( inputImage, 
    inputImage->GetLargestPossibleRegion() );  

  // the mask isn't required, so the iterator might not exist and need to be
  // iterated over
  MaskVolumePointer mask = this->GetMaskVolume();
  itk::ImageRegionConstIterator< MaskVolumeType > maskIt;

  bool isMaskUsed = false;
  if( mask ) {
    isMaskUsed = true;
  }

  if( isMaskUsed ) {
    maskIt = itk::ImageRegionConstIterator< MaskVolumeType >
      ( mask, mask->GetLargestPossibleRegion() );
  }
    
  if( isMaskUsed ) {
    maskIt = maskIt.Begin();
  }

  const int nTensorRows = 3;
  const int nTensorColumns = 3;

  // this initialization of the for loop is long
  for ( 
    inputImageIt = inputImageIt.Begin(),
    odfLookUpTableIt = odfLookUpTableIt.Begin();
    
    !inputImageIt.IsAtEnd() &&
    !odfLookUpTableIt.IsAtEnd();
    
    ++inputImageIt,
    ++odfLookUpTableIt )
  {
            
    bool isPixelMasked = false;
        
    if( isMaskUsed ) {
      isPixelMasked = maskIt.Value() == 0;
    }
        
    bool hasZero = false;        
    if( !isPixelMasked ) {
        
      // determine if a zero exists within the current tensor
      PixelType currentTensor = inputImageIt.Value();
      
      typedef typename PixelType::ConstIterator TensorIterator;
      for( TensorIterator tensorIterator = currentTensor.Begin(); 
        tensorIterator != currentTensor.End() && !hasZero; ++tensorIterator ) {
        
        if( *tensorIterator == 0.0 ) {
          hasZero = true;
        }
      
      }
          
    }
                        
    if( !hasZero && !isPixelMasked ) {

      /* DJ: my attempt at explaining what look is.  I'd like to get rid of the 
         look up table at some point:
         
         look is actually a 3D matrix, the size of the mask.  At every location
         that the volume contains a non-masked pixel, look will contain an index
         that cooresponds to the index that the odf is created in m_Odfs        

         MATLAB: look = zeros(size(mask));
      */

      odfLookUpTableIt.Set( cOdfs );          
      PixelType currentTensor = inputImageIt.Value();

      /* diffusion tensors stored with taking advantage of the storage savings
         offered by the symmmetry because we need to multiply the elements of
         the tensor by the polarity */       
      itk::Matrix< double, nTensorRows, nTensorColumns > tensor;

      // MATLAB: d = d(order, order) .* polarity;
      for( int cTensorRow=0; cTensorRow<nTensorRows; cTensorRow++ ) {
      
        // we want flip the tensors rows and columns backwards (not 
        // transpose though)
        const int cFlippedTensorRow = ( nTensorRows - 1 ) - cTensorRow;

        for( int cTensorColumn=0; cTensorColumn<nTensorColumns; cTensorColumn++ ) {

          const int cTensorFlippedColumn = ( nTensorColumns - 1 ) - cTensorColumn;
  
          const double flippedValue = currentTensor( cFlippedTensorRow, cTensorFlippedColumn );

// TODO: remove this if the energy is wrong...
//          const double flippedValue = currentTensor( cTensorRow, cTensorColumn );

          const double currentPolarity = polarity[ cTensorRow ][ cTensorColumn ];
                    
          tensor[ cTensorRow ][ cTensorColumn ] = flippedValue * currentPolarity;
//          tensor[ cTensorRow ][ cTensorColumn ] = currentTensor( cTensorRow, cTensorColumn );
          
        }
      
      }
  
      // MATLAB: d = R*d*invR;
      if( !isRotationIdentity ) {
        tensor = rotationMatrix * tensor.GetVnlMatrix() * inverseRotationMatrix;
      }
  
      ArrayPointer odf = new ArrayType( this->GetNumberOfDirections() );
      CalculateTensor2Odf( &tensor, odf);
                    
      this->m_Odfs.push_back( odf );
      
      if( ( cOdfs % 10000 ) == 0 ) {
        this->InvokeEvent( PoistatsOdfCalculationProgressEvent() );            
      }
      
      cOdfs++;
                
    }
    
    // increment the mask iterator if a mask is provided
    if( isMaskUsed ) {
      ++maskIt;    
    }
    
  }
  
  // stop the clock and calculat the elapsed time
  clock_t endClock = clock();
  const double elapsedTime = 
    static_cast< double >( endClock - startClock ) / CLOCKS_PER_SEC;
  this->SetElapsedTime( elapsedTime );
  this->InvokeEvent( PoistatsOdfCalculationEndEvent() );            
}

template <class TInputImage, class TOutputImage>
void
PoistatsFilter<TInputImage, TOutputImage>
::GetMagnetToSliceFrameRotation( MatrixPointer rotation ) {

  /* MATLAB:  
    sliceup = [normal_s normal_a 0];
    R = rot3u2v([1 0 0], sliceup);
  */

  const double identityValues[] = { 1, 0, 0 };  
  ArrayType identity( identityValues, 3 );
  
  if( this->GetInput() ) {

    const double normalS = this->GetInput()->GetDirection()( 2, 2 );
    const double normalA = this->GetInput()->GetDirection()( 1, 2 );
    const double normalR = 0.0;
    ArrayType sliceUp( 3 );
    sliceUp( 0 ) = normalS;
    sliceUp( 1 ) = normalA;
    sliceUp( 2 ) = normalR;
    
    GenerateRotationMatrix3u2v( &identity, &sliceUp, rotation );
    
  }
}

template <class TInputImage, class TOutputImage>
void
PoistatsFilter<TInputImage, TOutputImage>
::GenerateRotationMatrix3u2v( ArrayPointer u, ArrayPointer v, 
  MatrixPointer outputRotationMatrix ) {

  typedef vnl_matrix< double > VnlMatrixType;
  typedef vnl_vector< double > VnlVectorType;  

  //MATLAB: dotprod = u'*v;
  VnlMatrixType uMatrix( u->size(), 1 );
  uMatrix.set_column( 0, *u );

  VnlMatrixType vMatrix( v->size(), 1 );
  vMatrix.set_column( 0, *v );  
  
  VnlMatrixType dotProductMatrix = uMatrix.transpose() * vMatrix;
  const double dotProduct = dotProductMatrix[ 0 ][ 0 ];
  
// MATLAB:  
//  if abs(dotprod - 1) < 1e-3
//    % if vectors aligned
//    R = eye(3);    
//  else
  MatrixType rotation( u->size(), v->size() );

  const double tiny = 1e-3;
  if( fabs( dotProduct - 1 ) < tiny ) {
    // if vectors are aligned
    rotation.set_identity();
  } else {
// MATLAB:
//   % find rotation axis w
//  if abs(dotprod + 1) < 1e-3
//    % if vectors are anti-aligned
//    % rotation axis is undefined and
//    % so need to find arbitrary orthgonal axis
//    if u(1) == 0 & u(2) == 0
//      v2 = [1 0 0]';
//    else
//      v2 = [u(2) -u(1) 0]';  
//    end
//    v2 = v2/norm(v2);
//    w = cross(u,v2); 
//  else
//    w = cross(u,v); 
//  end

    // find rotation axis w
    if( (dotProduct + 1 ) < tiny ) {
      // if vectors are anti-aligned, rotation axis is undefined and so need 
      //  to find arbitrary orthogonal axis
      
      if( ( uMatrix[ 0 ][ 0 ] == 0 ) && ( uMatrix[ 1 ][ 0 ] == 0 ) ) {
        const double identity[] = { 1, 0, 0 };
        vMatrix.set_column( 0, identity );
      } else {
        const double newV[] = { uMatrix[ 1 ][ 0 ], -uMatrix[ 0 ][ 0 ], 0.0 };
        vMatrix.set_column( 0, newV );
      }
      
      const double twoNorm = vMatrix.array_two_norm();
      vMatrix /= twoNorm;
    }
      
    VnlVectorType uVector = uMatrix.get_column( 0 );
    VnlVectorType vVector = vMatrix.get_column( 0 );

    VnlVectorType wVector = vnl_cross_3d< double >( uVector, vVector );

    // MATLAB:    
    //    w = w/norm(w);
    //    x = w(1); 
    //    y = w(2); 
    //    z = w(3);
    VnlMatrixType wMatrix( 3, 1 );
    wMatrix.set_column( 0, wVector );
    const double twoNorm = wMatrix.array_two_norm();
    wMatrix /= twoNorm;
    
    const double x = wMatrix[ 0 ][ 0 ];
    const double y = wMatrix[ 1 ][ 0 ];
    const double z = wMatrix[ 2 ][ 0 ];

    // MATLAB:  
    //    P = [+0 -z +y;
    //         +z +0 -x;
    //         -y +x +0];         

    const double pValues[] = {
       0, -z,  y,
       z,  0, -x,
      -y,  x,  0
    };
    VnlMatrixType p( pValues, 3, 3 );

    // MATLAB:
    //    theta = acos(u'*v);
    VnlMatrixType angle = uMatrix.transpose() * vMatrix;
    const double theta = acos( angle[ 0 ][ 0 ] );

    // MATLAB: R = eye(3) + P*sin(theta) + P*P*(1-cos(theta));

    VnlMatrixType identity( 3, 3 );
    identity.set_identity();
    
    rotation = identity + p * sin( theta ) + p * p * ( 1 - cos (theta ) );
  }

  // copy the matrix to the output  
  for( unsigned int cRow=0; cRow<outputRotationMatrix->rows(); cRow++ ) {
    for( unsigned int cCol=0; cCol<outputRotationMatrix->cols(); cCol++ ) {
      ( *outputRotationMatrix )[ cRow ][ cCol ] = rotation[ cRow ][ cCol ];
    }
  }
  
}

template <class TInputImage, class TOutputImage>
void
PoistatsFilter<TInputImage, TOutputImage>
::RoundPath( itk::Array2D< int > *outputRoundedArray, 
             itk::Array2D< double >* inputArray ) {
  
  for( unsigned int cRow = 0; cRow < inputArray->rows(); cRow++ ) {
    for( unsigned int cColumn = 0; cColumn < inputArray->cols(); cColumn++ ) {
    
      ( *outputRoundedArray )[ cRow ][ cColumn] = 
        static_cast< int >( round( ( *inputArray)[ cRow ][ cColumn] ) );
      
    }
  }
  
}

template <class TInputImage, class TOutputImage>
typename PoistatsFilter<TInputImage, TOutputImage>::OutputImageType *
PoistatsFilter<TInputImage, TOutputImage>
::GetOptimalDensity() {
  return dynamic_cast<OutputImageType*>( 
    this->ProcessObject::GetOutput(PATH_DENSITY_OUTPUT) );
}


template <class TInputImage, class TOutputImage>
typename PoistatsFilter<TInputImage, TOutputImage>::OutputImageType *
PoistatsFilter<TInputImage, TOutputImage>
::GetDensity() {
  return dynamic_cast<OutputImageType*>( 
    this->ProcessObject::GetOutput(OPTIMAL_PATH_DENSITY_OUTPUT) );
}

template <class TInputImage, class TOutputImage>
int
PoistatsFilter<TInputImage, TOutputImage>
::GetNumberOfInitialPoints() const {
  return this->m_Seeds.size();
}

/**
 * Returns the path with evenly spaced sample points.
 */
template <class TInputImage, class TOutputImage>
void
PoistatsFilter<TInputImage, TOutputImage>
::InitPaths() {
  /* MATLAB:
    % initialize paths
    fprintf('Initializing paths ...\n');  
    origpath = rethreadpath(initpoints, ncontrolpoints+2);
    lowtrialpath = origpath;
    for i = 1:nreplica
      basepath{i} = origpath;
      prevpath{i} = rethreadpath(origpath, steps);
      trialpath{i} = zeros(steps,3);
      bestpath{i} = zeros(steps,3);
    end  
  */

  if( this->m_PoistatsModel->IsUsingPathInitialization() ) {
    this->m_Replicas->InitializePathsUsingEigenVectors();
  } else {

    const int nInitialPoints = this->GetNumberOfInitialPoints();
    const int spatialDimensions = 3;
    
    MatrixType initialPoints( nInitialPoints, spatialDimensions );
    this->GetInitialPoints( &initialPoints );
  
    this->m_Replicas->SetInitialPoints( &initialPoints );
  }
  
  this->InvokeEvent( SeedsFoundInitialEvent() );
}

/**
 * Returns the starting end points based on the seed locations.
 */
template <class TInputImage, class TOutputImage>
void
PoistatsFilter<TInputImage, TOutputImage>
::GetInitialPoints( MatrixPointer initialPoints ) {

  /* MATLAB:
    for idx = 1:nseeds
      [i j k] = ind2sub(size(seeds),find(seeds==seedvalues(idx)));
      X = [i j k];
      com = mean(X,1); % com = center of mass
      % jdx = subscript into voxel distance list of smallest distance
      % between com and each voxel in region
      [null jdx] = min(distmat(com, X));
      % min can return multiples in case of ties. arbitrarily select first tie.
      jdx = jdx(1);
      % initpoints = dim1 x dim2 array of subscripts of voxels closest to coms
      % of each roi. dim1 = number of seeds. dim2 = 3 (image dimension)
      initpoints(idx,:) = X(jdx,:);  
    end
  */

  const int nDimensions = 3;

  ArrayType closestPoint( nDimensions );
  for( unsigned int cSeed=0; cSeed<this->m_Seeds.size(); cSeed++ ) {

    MatrixPointer seedRegion = this->m_Seeds[ cSeed ];
    
    GetPointClosestToCenter( seedRegion, &closestPoint );    
    itkDebugMacro( << "closest point: " << closestPoint << std::endl );
    
    if( cSeed < initialPoints->rows() ) {
    
      for( int cDim=0; cDim<nDimensions; cDim++ ) {
        ( *initialPoints )[ cSeed ][ cDim ] = closestPoint[ cDim ];
      }
      
    } else {
    
      itkGenericExceptionMacro (<< "too many seeds");            
      
    }
    
  }
  
}

template <class TInputImage, class TOutputImage>
void
PoistatsFilter<TInputImage, TOutputImage>
::GetPointClosestToCenter( MatrixPointer seeds, ArrayPointer closestSeed ) {

  const int nDimensions=closestSeed->size();
  
  ArrayType center( nDimensions );
  GetCenterOfMass( seeds, &center );
  
  double minDistance = std::numeric_limits< double >::max();
  int indexOfMin = 0;
  
  for( unsigned int cSeed=0; cSeed<seeds->rows(); cSeed++ ) {  
    vnl_vector< double > seed = seeds->get_row( cSeed );
    
    double distance = GetDistance( &seed, &center );

    // if this distance is the least, then save it
    if( distance < minDistance ) {
      minDistance = distance;
      indexOfMin = cSeed;
    }
    
  }
  
  // copy the closest seed over
  for( int cDim=0; cDim<nDimensions; cDim++ ) {
    ( *closestSeed )[ cDim ] = ( *seeds )[ indexOfMin ][ cDim ];
  }

}

template <class TInputImage, class TOutputImage>
double
PoistatsFilter<TInputImage, TOutputImage>
::GetDistance( vnl_vector< double >* point1, ArrayPointer point2 ) {

  // subtract them element by element
  vnl_vector< double > distanceVector = *point1 - *point2;
  
  double distance = 0.0;
  
  // square each element
  for( unsigned int cDim=0; cDim<distanceVector.size(); cDim++ ) {
    distance += distanceVector[ cDim ] * distanceVector[ cDim ];
  }
  
  distance = sqrt( distance );
  return distance;
}

template <class TInputImage, class TOutputImage>
void
PoistatsFilter<TInputImage, TOutputImage>
::GetCenterOfMass( MatrixPointer mass, ArrayPointer center ) {

  const double nSeeds = static_cast< double >( mass->rows() );
  const int nDimensions = mass->cols();
  
  // get the mean along the columns
  for( int cDim=0; cDim<nDimensions; cDim++ ) {

    // initialize
    ( *center )[ cDim ] = 0.0;

    for( unsigned int cSeeds=0; cSeeds<mass->rows(); cSeeds++ ) {
        
      ( *center )[ cDim ] += ( *mass )[ cSeeds ][ cDim ];
    
    }
    
    ( *center )[ cDim ] /= nSeeds;
    
  }
  
}

template <class TInputImage, class TOutputImage>
int
PoistatsFilter<TInputImage, TOutputImage>
::GetNumberOfControlPoints() {
  return this->m_Replicas->GetNumberOfControlPoints();
}

template <class TInputImage, class TOutputImage>
void
PoistatsFilter<TInputImage, TOutputImage>
::SetNumberOfControlPoints( const int nPoints ) {
  this->m_Replicas->SetNumberOfControlPoints( nPoints );
}

template <class TInputImage, class TOutputImage>
void
PoistatsFilter<TInputImage, TOutputImage>
::SetNumberOfSteps( const int nSteps ) {
  this->m_Replicas->SetNumberOfSteps( nSteps );
}

template <class TInputImage, class TOutputImage>
int 
PoistatsFilter<TInputImage, TOutputImage>
::GetNumberOfSteps() {
  return this->m_Replicas->GetNumberOfSteps();
}

template <class TInputImage, class TOutputImage>
void
PoistatsFilter<TInputImage, TOutputImage>
::CalculateTensor2Odf( itk::Matrix< double, 3, 3 > *tensor,
                       itk::Array< double > *odf ) {

  const int numberOfOdfs = this->GetNumberOfDirections();
  const int numberOfTensorAxes = 3;  

  // this is the matlab operations that we're doing
  //    odf = sum((geo*pinv(D)).*geo,2).^(-3/2);

  VnlMatrixType geo = this->GetTensorGeometry();

  vnl_svd< double > tensorSvd( tensor->GetVnlMatrix() );
  VnlMatrixType tensorInverse( tensorSvd.pinverse( 3 ) );

  VnlMatrixType morphedGeo( geo * tensorInverse );
  
  typedef vcl_complex< double > ComplexType;
  typedef itk::Array< ComplexType > ComplexArray;
  ComplexArray complexOdf( odf->size() );
  
  ComplexType odfSum( 0.0, 0.0 );
  for( int cRow=0; cRow<numberOfOdfs; cRow++ ) {
    
    // reset the sum for every odf
    double rowSum = 0.0;
    
    // all the axes (x,y,z) of this odf direction for the sum of this direction
    for( int cColumn=0; cColumn<numberOfTensorAxes; cColumn++ ) {
      rowSum += morphedGeo[ cRow ][ cColumn ] * geo[ cRow ][ cColumn ];
    }

    // save the sum of this particular direction
        
    ComplexType complexRowSum( rowSum, 0 );
    complexOdf[ cRow ] = vcl_pow( complexRowSum, -1.5 );
    
    odfSum += complexOdf[ cRow ];

  }
  
  // normalize the odf
  //  odf = odf/sum(odf);
  for( int cOdf=0; cOdf<numberOfOdfs; cOdf++ ) {    
    ( *odf )[ cOdf ] = vcl_abs( complexOdf[ cOdf ] / odfSum );
  }
  
}

template <class TInputImage, class TOutputImage>
double
PoistatsFilter<TInputImage, TOutputImage>
::CalculateOdfPathEnergy( MatrixPointer path,
                          ArrayPointer* odfs,
                          ArrayPointer outputEnergies ) {
  
  const int spatialDimension = path->columns();
  
  MatrixType pathDifference( path->rows()-1, spatialDimension );
  PoistatsReplica::CalculatePathVectors( path, &pathDifference );
  
  ArrayType magnitude( path->rows() - 1 );
  PoistatsReplica::CalculateMagnitude( &pathDifference, &magnitude );  

  MatrixType normalizedPathVectors( pathDifference );

  // normalize the vectors
  for( unsigned int cPath=0; cPath<pathDifference.rows(); cPath++ ) {
    for( unsigned int cDimension=0; cDimension<pathDifference.cols(); cDimension++ ) {
      normalizedPathVectors[ cPath ][ cDimension ] /= magnitude[ cPath ];
    }
  }

  ArrayType anglesBetweenPathVectors( normalizedPathVectors.rows()-1 );
  CalculateAnglesBetweenVectors( &normalizedPathVectors, 
                                 &anglesBetweenPathVectors );

  // MATLAB: if max(angles) > pi/3, energies = 1e6; end;
  static const double largeAngle = M_PI / 3.0;    
  double meanEnergy = 0.0;
  
  // we don't want sharp turns, so set those very large
  if( anglesBetweenPathVectors.max_value() > largeAngle ) {
  
    const double maxEnergy = 1e6;
    meanEnergy = maxEnergy / magnitude.sum();

    if( outputEnergies ) {
      outputEnergies->Fill( maxEnergy );
    }

  } else {
    
    // find dot products between tangent vectors and geometry points    
    vnl_matrix< double > geo = this->GetTensorGeometry();
    
    // in order to get the exact same results as the matlab version, I need to
    // swap the first and last columns
    for( unsigned int row=0; row<geo.rows(); row++ ) {
      const double tmp = geo[ row ][ 0 ];
      geo[ row ][ 0 ] = geo[ row ][ 2 ];
      geo[ row ][ 2 ] = tmp;
    }
    
    vnl_matrix< double > dotProductPerGeoDirection( normalizedPathVectors * 
                                                    geo.transpose() );

    // take the absolute values of the dot products
    for( unsigned int cRow=0; cRow<dotProductPerGeoDirection.rows(); cRow++ ) {
      for( unsigned int cColumn=0; cColumn<dotProductPerGeoDirection.cols(); cColumn++ ) {
        dotProductPerGeoDirection[ cRow ][ cColumn ] = 
          fabs( dotProductPerGeoDirection[ cRow ][ cColumn ] );
      }
    }
  
    // calculate the angles between path and geo
    MatrixType pathGeoAngles( dotProductPerGeoDirection.rows(), dotProductPerGeoDirection.cols() );
    for( unsigned int cRow=0; cRow<pathGeoAngles.rows(); cRow++ ) {
      for( unsigned int cColumn=0; cColumn<pathGeoAngles.cols(); cColumn++ ) {
        pathGeoAngles[ cRow ][ cColumn ] = 
          acos( dotProductPerGeoDirection[ cRow ][ cColumn ] );
      }
    }  
        
    MatrixType densityMatrix( pathGeoAngles.rows(), pathGeoAngles.cols() );
    CalculateDensityMatrix( &pathGeoAngles, &densityMatrix );
        
    /* MATLAB:
      odflist = odflist(1:(end-1),:); % n - 1 sets of odfs, because they
                                    % correspond to path segments, not path points.    
    % odfvalues are the baysian conditional posterior probability distribution
    %  sum(A,2) will sum along the rows
    odfvalues = sum(odflist.*dm,2);
    */
    ArrayType odfValues( densityMatrix.rows() );
    for( unsigned int cRow=0; cRow<densityMatrix.rows(); cRow++ ) {
      
      double odfListSum = 0.0;
      
      // sums along the rows
      for( unsigned int cColumn=0; cColumn<densityMatrix.cols(); cColumn++ ) {
        ArrayPointer odfsAtRow = odfs[ cRow ];
        const double odf = ( *odfsAtRow )[ cColumn ] ;
        const double density = densityMatrix[ cRow ][ cColumn ];
        odfListSum += odf * density;
      }
      
      odfValues[ cRow ] = odfListSum;
    }
    
    /* MATLAB:    
    % compute energy
    % the energy is defined as the negative logarithm of the conditional
    % confirmational posterior distribution see Habeck et el Physical Review E
    % 72, 031912, 2005, equation 17
    energies = -log(abs(odfvalues)).*nm;
    */
    ArrayType energies( odfValues.size() );
    for( unsigned int cRow=0; cRow<energies.size(); cRow++ ) {
      energies[ cRow ] = -log( fabs( odfValues[ cRow ] ) ) * magnitude[ cRow ];
    }
    
    // if the output energies exist, we should copy them out
    if( outputEnergies ) {
      for( unsigned int cRow=0; cRow<energies.size(); cRow++ ) {
        ( *outputEnergies )[ cRow ] = energies[ cRow ];
      }
    }
  
    // MATLAB: meanenergy = sum(energies)/sum(nm);
    meanEnergy = energies.sum() / magnitude.sum();
  }
   
  return meanEnergy;
}

template <class TInputImage, class TOutputImage>
void
PoistatsFilter<TInputImage, TOutputImage>
::CalculateDensityMatrix( itk::Array2D< double > *angles, 
  itk::Array2D< double > *densityMatrix ) {

// dm = rownormalize(exp(-8*angles.^2));

  for( unsigned int cRow=0; cRow<angles->rows(); cRow++ ) {

    double rowSum = 0.0;

// this was the default
//    const double variance = 8;

// this worked well when I initialized all the points on the helix
//    const double variance = 8*8;

    const double variance = 8*6;

    for( unsigned int cColumn=0; cColumn<angles->cols(); cColumn++ ) {

      const double angle = ( *angles )[ cRow ][ cColumn ];
      const double currentDensity = exp( -variance * angle * angle );
      
      ( *densityMatrix )[ cRow ][ cColumn ] = currentDensity;
      rowSum += currentDensity;

    }
    
    // normalize based on the row
    for( unsigned int cColumn=0; cColumn<angles->cols(); cColumn++ ) {
      ( *densityMatrix )[ cRow ][ cColumn ] /= rowSum;
    }
    
  }
}

//% angles are the angles between adjacent path vectors, v
//%   ie. v1 dot v2 = cos(theta12)
//angles = acos(abs(sum((v.*circshift(v,[-1 0])),2)));
template <class TInputImage, class TOutputImage>
void
PoistatsFilter<TInputImage, TOutputImage>
::CalculateAnglesBetweenVectors( itk::Array2D< double > *vectors, 
  itk::Array< double > *angles ) {

  for( unsigned int cRow=0; cRow<vectors->rows()-1; cRow++ ) {

    double rowSum = 0.0;

    for( unsigned int cColumn=0; cColumn<vectors->cols(); cColumn++ ) {

      const double currentCell = ( *vectors )[ cRow ][ cColumn ];
      const double nextCell = ( *vectors )[ cRow+1 ][ cColumn ];
      rowSum += currentCell * nextCell;

    }

    ( *angles )[ cRow ] = acos( fabs( rowSum ) );

  }  

}

template <class TInputImage, class TOutputImage>
void 
PoistatsFilter<TInputImage, TOutputImage>
::GetOdfsAtPoints( itk::Array< double >** outputOdfs, 
                   itk::Array2D< int >* inputPoints ) {
//  idx = inboundsidx(rpath, size(mask));
      
//  % now get the indices that aren't masked
//  odfidx = look(idx);
//  goodindices = odfidx~=0;     
//  odfs(goodindices,:) = abs(odflist(odfidx(goodindices),:));      

  InputImageConstPointer inputImage = this->GetInput();

  SizeType imageSize = this->GetInput()->GetLargestPossibleRegion().GetSize();

  for( unsigned int cPoint=0; cPoint<inputPoints->rows(); cPoint++ ) {

    OdfLookUpIndexType index;

    index[ 0 ] = ( *inputPoints )[ cPoint ][ 0 ];
    index[ 1 ] = ( *inputPoints )[ cPoint ][ 1 ];
    index[ 2 ] = ( *inputPoints )[ cPoint ][ 2 ];
    
    bool isValidIndex = true;                
    for( unsigned int cDim=0; cDim<OdfLookUpRegionType::GetImageDimension(); cDim++ ) {
      if( index[ cDim ] < 0 ) {
        isValidIndex = false;
      } else if( static_cast< unsigned int >( index[ cDim ] ) >= imageSize[ cDim ] ) {
        isValidIndex = false;
      }
    }
    
    OdfLookUpTablePointer table = this->GetOdfLookUpTable();
    int odfIndex = INVALID_INDEX;
    if( isValidIndex ) {
      odfIndex = table->GetPixel( index );
    }
    
    if( odfIndex > INVALID_INDEX ) {    
      outputOdfs[ cPoint ] = this->m_Odfs[ odfIndex ];            
    } else {
      outputOdfs[ cPoint ] = &m_InvalidOdf;
    }
  }

}

template <class TInputImage, class TOutputImage>
void 
PoistatsFilter<TInputImage, TOutputImage>
::GetSortedUniqueSeedValues( SeedVolumePointer volume, 
  std::vector< std::pair< SeedType, int > > *seedValues ) {

  typedef std::pair< SeedType, int > SeedValueCountPairType;
  typedef std::vector< SeedValueCountPairType > SeedPairList;

  itk::ImageRegionConstIterator< SeedVolumeType > seedIt(
    volume, volume->GetLargestPossibleRegion() );

  // get the unique values
  for ( seedIt = seedIt.Begin(); !seedIt.IsAtEnd(); ++seedIt ) {
    
    const SeedType pixelValue = seedIt.Value();
    
    if( pixelValue != 0 ) {
 
      bool isUnique = true;
      
      // where the seed value should be inserted to maintain sorted order
      SeedPairList::iterator insertionIndex = seedValues->begin();
  
      for( SeedPairList::iterator valuesIt = seedValues->begin();
          valuesIt != seedValues->end(); valuesIt++ ) {
        
        const SeedType seedValue = ( *valuesIt ).first;
        
        if( pixelValue == seedValue ) {
          isUnique = false;
          // increment the count
          ( *valuesIt ).second++; 
        } else if( pixelValue > seedValue ) {
          insertionIndex = valuesIt;
        }
        
      }
      
      if( isUnique ) {
        const int size = 1;
        SeedValueCountPairType seedPair( pixelValue, size );
        seedValues->insert( insertionIndex, seedPair );
      }

    }
  }
  
}

template <class TInputImage, class TOutputImage>
itk::Array2D< double >*
PoistatsFilter<TInputImage, TOutputImage>
::GetStartSeeds() {

  return this->m_Seeds[ 0 ];
}

template <class TInputImage, class TOutputImage>
itk::Array2D< double >*
PoistatsFilter<TInputImage, TOutputImage>
::GetEndSeeds() {
  return this->m_Seeds[ this->m_Seeds.size()-1 ];
}

template <class TInputImage, class TOutputImage>
itk::Array< double >
PoistatsFilter<TInputImage, TOutputImage>
::GetSamples() {
  
  typedef itk::BSplineInterpolateImageFunction< SamplingVolumeType, double, double > InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage( this->m_SamplingVolume );
  interpolator->SetSplineOrder( 3 );

  MatrixType finalPath = this->GetFinalPath();
        
  ArrayType samples( finalPath.rows() );  
  typedef InterpolatorType::ContinuousIndexType ContinuousIndexType;  
  ContinuousIndexType index;
  // evaluate at the final path points
  for( unsigned int cPoint=0; cPoint<finalPath.rows(); cPoint++ ) {
    
    for( unsigned int cDim=0; cDim<finalPath.cols(); cDim++ ) {
      index[ cDim ] = finalPath[ cPoint ][ cDim ];
    }
    
    samples[ cPoint ] = interpolator->EvaluateAtContinuousIndex( index );
  }
      
  return samples;
}

/**
 * Return in seeds1 the union of seeds1 and seeds2.  Throws an exception if
 * all seeds in seeds2 aren't found.
 */
template <class TInputImage, class TOutputImage>
void 
PoistatsFilter<TInputImage, TOutputImage>
::TakeUnionOfSeeds( std::vector< std::pair< SeedType, int > > *seeds1,
                    std::vector< SeedType > *seeds2 ) {
                    
  typedef std::pair< SeedType, int > SeedValueCountPairType;
  typedef std::vector< SeedValueCountPairType > SeedPairList;
  
  SeedPairList seedUnion;
  
  for( std::vector< SeedType >::iterator values2It = seeds2->begin();
    values2It != seeds2->end(); ++values2It ) {
        
    const SeedType value2 = ( *values2It );
    
    bool isSeedFound = false;
    
    for( SeedPairList::iterator values1It = seeds1->begin();
      values1It != seeds1->end(); ++values1It ) {
  
      const SeedType value1 = ( *values1It ).first;

      if( value1 == value2 ) {
        seedUnion.push_back( *values1It );
        isSeedFound = true;
      }

    }
    
    if( !isSeedFound ) {
      itkGenericExceptionMacro (<< "seed " << value2 << " not found in seed volume.");    
    }
            
  }
  
  seeds1->clear();
  *seeds1 = seedUnion;

}


/**
 * Set the seed volume and parses the volume for the number of different seeds
 * and the pixel values for those seeds.
 */
template <class TInputImage, class TOutputImage>
void 
PoistatsFilter<TInputImage, TOutputImage>
::SetSeedVolume( SeedVolumePointer volume ) {

  this->m_SeedVolume = volume;

}


template <class TInputImage, class TOutputImage>
void 
PoistatsFilter<TInputImage, TOutputImage>
::ParseSeedVolume() {
  
  typedef std::pair< SeedType, int > SeedValueCountPairType;
  typedef std::vector< SeedValueCountPairType > SeedPairList;
  SeedPairList seedValueCountPairs;
  GetSortedUniqueSeedValues( this->m_SeedVolume, &seedValueCountPairs );
  
  if( this->m_SeedValuesToUse.empty() ) {
    this->InvokeEvent( SeedsUsingAllEvent() );
  } else {
    try {
      TakeUnionOfSeeds( &seedValueCountPairs, &this->m_SeedValuesToUse );
    } catch( itk::ExceptionObject & excp ) {
      throw excp;
    }
  }
  
  itk::ImageRegionConstIterator< SeedVolumeType > seedIt(
    this->m_SeedVolume, this->m_SeedVolume->GetLargestPossibleRegion() );

  if( seedValueCountPairs.size() > 1 ) {
  
    // save the seed values that will be used in the model
    std::vector< int > *seedValues = new std::vector< int >;
    for( SeedPairList::iterator valuesIt = seedValueCountPairs.begin();
      valuesIt != seedValueCountPairs.end(); valuesIt++ ) {
                
      const SeedType seedValue = ( *valuesIt ).first;
      seedValues->push_back( seedValue );
      
    }
    this->m_PoistatsModel->SetSeedValues( seedValues );
  
    // save the pixel values
    for( SeedPairList::iterator valuesIt = seedValueCountPairs.begin();
      valuesIt != seedValueCountPairs.end(); valuesIt++ ) {
          
      const SeedType seedValue = ( *valuesIt ).first;
      const int nCurrentSeed = ( *valuesIt ).second;
      
      itkDebugMacro( << "  ( seed value, number of seeds ): (" << seedValue << ", " << nCurrentSeed << " )" );
  
      MatrixPointer currentSeeds = 
        new MatrixType( nCurrentSeed, SeedVolumeIndexType::GetIndexDimension() );
  
      int cCurrentSeeds = 0;
  
      this->m_Seeds.push_back( currentSeeds );
  
      for ( seedIt = seedIt.Begin(); !seedIt.IsAtEnd(); ++seedIt ) {    
  
        SeedType pixelValue = seedIt.Value();
  
        if( pixelValue == seedValue ) {
 
          // this is the seed index that we save
          SeedVolumeIndexType seedIndex = seedIt.GetIndex();
          
          // this is a long for loop
          for( unsigned int cIndex=0; cIndex<SeedVolumeIndexType::GetIndexDimension(); 
            cIndex++ ) {
            
            ( *currentSeeds )[ cCurrentSeeds  ][ cIndex  ] = seedIndex[ cIndex ];
            
          }
  
          cCurrentSeeds++;
          
        }
        
      }
      
    }
  
  } else {
    itkGenericExceptionMacro (<< "not enough seed regions.");    
  }

}

template <class TInputImage, class TOutputImage>
void
PoistatsFilter<TInputImage, TOutputImage>
::CalculateFinalPathProbabilities() {

  ArrayType probabilities( this->GetNumberOfSamplePoints() );

  MatrixType finalPath = this->GetFinalPath();

  this->GetPathProbabilities( &finalPath, &probabilities );
  
  this->m_FinalPathProbabilities = probabilities;

}

template <class TInputImage, class TOutputImage>
void
PoistatsFilter<TInputImage, TOutputImage>
::GetPathProbabilities( MatrixPointer path, ArrayPointer probabilities ) {

/*
odfs = (1e-200)*ones(steps, length(geo));
rpath = round(globalbestpath);      
idx = inboundsidx(rpath, size(mask));
odfidx = look(idx);
goodindices = odfidx~=0;     
odfs(goodindices,:) = abs(odflist(odfidx(goodindices),:));
[energy energies] = odfpathenergy(globalbestpath, odfs, geo);             
pathprobabilities = exp(-energies);
pathprobabilities = csapi(1:length(pathprobabilities), pathprobabilities, ...
                          linspace(1,length(pathprobabilities),nsamplepoints));
*/

  ArrayPointer odfs[ path->rows() ];
                   
  itk::Array2D< int > roundedPath( path->rows(), path->cols() );                                 

  RoundPath( &roundedPath, path );

  this->GetOdfsAtPoints( odfs, &roundedPath );
  
  ArrayType energies( path->rows()-1 );

  CalculateOdfPathEnergy( path, odfs, &energies );
  
  MatrixType rawProbablities( energies.size(), 3 );
  rawProbablities.Fill( 0.0 );

  for( unsigned int cProbability=0; cProbability<rawProbablities.rows(); 
    cProbability++ ) {
    
    rawProbablities[ cProbability ][ 0 ] = exp( -energies[ cProbability ] );
      
  }
  
  const double gridFloor = 0.0;
  const double gridCeiling = 1.0;
  ArrayType originalPathGrid( rawProbablities.rows() );
  PoistatsReplica::SpaceEvenly( &originalPathGrid, gridFloor, gridCeiling );  
  
  MatrixPointer outputProbabilities = m_Replicas->CubicSplineInterpolation( 
    &rawProbablities,
    &originalPathGrid,
    probabilities->size() );
    
  for( unsigned int cProbability=0; cProbability<probabilities->size(); 
    cProbability++ ) {
    
    ( *probabilities )[ cProbability ] = 
      ( *outputProbabilities )[ cProbability ][ 0 ];
  
  }
  
}

template <class TInputImage, class TOutputImage>
void
PoistatsFilter<TInputImage, TOutputImage>
::SetRandomSeed( const long seed ) {

  this->m_PoistatsModel->SetRandomSeed( seed );

}


template <class TInputImage, class TOutputImage>
void
PoistatsFilter<TInputImage, TOutputImage>
::GetAggregateReplicaDensities( MatrixListType replicaPaths, 
  OutputImagePointer aggregateDensities ) {

  /* MATLAB:
  ps = exp(-energy); ps = ps/sum(ps);
  density = zeros(size(mask));
  for i = 1:nreplica
    if ~mod(i,5), fprintf('%d ',i);end;
    if ps(i) > (1e-1/nreplica) % if replica contributes
      density = density + ps(i)*pnts2vol(bestpath{i}, size(density),.5);
    end
  end
  density = density / sum(density(:));
  */

  ArrayType probability( m_Replicas->GetNumberOfReplicas() );
  for( unsigned int cProbability=0; cProbability<probability.size(); cProbability++ ) {
    probability[ cProbability ] = 
      exp( -m_Replicas->GetCurrentMeanEnergy( cProbability ) );
  }  
  probability /= probability.sum();

  aggregateDensities->FillBuffer( 0.0 );
  
  OutputSizeType size;
  OutputIndexType start;

  for( unsigned int cDim=0; cDim<OutputRegionType::GetImageDimension(); cDim++ ) {

    size[ cDim ] = this->GetInput()->GetLargestPossibleRegion().GetSize( cDim );
    start[ cDim ] = 0;
  
  }
  
  OutputRegionType region;
  region.SetSize( size );
  region.SetIndex( start );
  
  typedef itk::ImageRegionIterator< OutputImageType > DensityIteratorType;

  DensityIteratorType aggregateIterator(
    aggregateDensities, aggregateDensities->GetLargestPossibleRegion() );
  
  const double minimumReplicaContributation = 
    1e-1 / static_cast< double >( replicaPaths.size() );
    
  int cReplica = 0;
  for( MatrixListType::iterator replicaIterator = replicaPaths.begin(); 
    replicaIterator != replicaPaths.end(); replicaIterator++ ) {

    this->InvokeEvent( GenerateBestReplicaPathDensitiesProgressEvent() );
    const MatrixPointer currentReplicaPath = ( *replicaIterator );
        
    OutputImagePointer currentReplicaDensity = OutputImageType::New();
    currentReplicaDensity->SetRegions( region );
    currentReplicaDensity->Allocate();
    currentReplicaDensity->FillBuffer( 0 );
    
    const double currentReplicaProbablity = probability[ cReplica ];
    
    if( currentReplicaProbablity > minimumReplicaContributation ) {
  
      this->ConvertPointsToImage( currentReplicaPath, currentReplicaDensity);
      DensityIteratorType currentReplicaIterator(
        currentReplicaDensity,
        currentReplicaDensity->GetLargestPossibleRegion() );

      for(
        currentReplicaIterator = currentReplicaIterator.Begin(),
        aggregateIterator = aggregateIterator.Begin();
        
        !currentReplicaIterator.IsAtEnd() &&
        !aggregateIterator.IsAtEnd();
        
        ++currentReplicaIterator,
        ++aggregateIterator ) {
        
        const double currentReplicaContribution = 
          currentReplicaIterator.Value() * currentReplicaProbablity;

        aggregateIterator.Set( aggregateIterator.Value() + 
          currentReplicaContribution );
        
      }
    
    }
    
    cReplica++;
    
  }

  // get the sum of the image
  double aggregateSum = 0.0;  
  for( aggregateIterator = aggregateIterator.Begin(); 
    !aggregateIterator.IsAtEnd(); ++aggregateIterator ) {
    aggregateSum += aggregateIterator.Value();
  }

  // normalize the image
  const double inverseAggregateSum = 1 / aggregateSum;
  for( aggregateIterator = aggregateIterator.Begin(); 
    !aggregateIterator.IsAtEnd(); ++aggregateIterator ) {
    aggregateIterator.Set( aggregateIterator.Value() * inverseAggregateSum );
  }
  
}

/**
 * Sets the maks volume to use.  This is optional.  If set, all non-zero areas
 * are masked when finding paths.  The areas that the path can go will be marked
 * with 1 and where the paths can't go will be 0.
 * 
 * @param volume Mask volume.
 */
template <class TInputImage, class TOutputImage>
void
PoistatsFilter<TInputImage, TOutputImage>
::SetMaskVolume( MaskVolumePointer volume ) {
  
  this->m_MaskVolume = volume;
  
  typedef itk::ImageRegionIterator< MaskVolumeType > MaskIteratorType;
  MaskIteratorType maskIterator(
    this->m_MaskVolume, this->m_MaskVolume->GetLargestPossibleRegion() );

  const MaskType validValue = 1;
  const MaskType invalidValue = 0;

  // the mask that the user inputs will have non zero areas for places that the
  // path can go.  In memory, we'd like 0 to be be places that the mask can't
  // go 1 to be where the path can go.
  
  // go through all pixels and assign value of 1 for all non-zero
  for( maskIterator = maskIterator.Begin(); !maskIterator.IsAtEnd(); ++maskIterator ) {
    
    MaskType value = maskIterator.Value();
        
    if( value != invalidValue && value != validValue ) {
      maskIterator.Set( validValue );
    }
    
  }

}

/**
 * Sets the next seed value to be used in connecting the path.
 * @param seedValue Next seed value to be used--pushed to the vector.
 */
template <class TInputImage, class TOutputImage>
void
PoistatsFilter<TInputImage, TOutputImage>
::SetNextSeedValueToUse( int seedValue ) {
  this->m_SeedValuesToUse.push_back( seedValue );
}

template <class TInputImage, class TOutputImage>
typename PoistatsFilter<TInputImage, TOutputImage>::VnlMatrixType
PoistatsFilter<TInputImage, TOutputImage>
::GetTensorGeometry() {

  if( this->m_TensorGeometry.empty() ) {
  
    const int numberOfOdfs = this->GetNumberOfDirections();
    const int numberOfTensorAxes = 3;  

    this->m_TensorGeometry = VnlMatrixType( *NO_ZERO_SHELL_252,
      numberOfOdfs, numberOfTensorAxes );
          
  }
  
  return this->m_TensorGeometry;
  
}

/**
 * Takes the union of the mask volume (if it exists) and the seed regions.  This
 * is done in case the mask doesn't include the seed volume.
 */
template <class TInputImage, class TOutputImage>
void
PoistatsFilter<TInputImage, TOutputImage>
::TakeUnionOfMaskAndSeeds() {

  bool isUnmasked = false;

  MaskVolumePointer mask = this->GetMaskVolume();
  
  // if the mask exists, take the union
  if( mask ) {
  
    // go through all the seed regions and make sure that it's a valid region
    // in the mask
    for( unsigned int cSeed=0; cSeed<this->m_Seeds.size(); cSeed++ ) {  

      MatrixPointer seedRegion = this->m_Seeds[ cSeed ];
      
      // go through each point in a region and make sure it's valid within the
      // mask
      for( unsigned int row=0; row<seedRegion->rows(); row++ ) {
      
        // go through each dimension of the seed and create a mask index
        vnl_vector< double > seedIndex = seedRegion->get_row( cSeed );        
        MaskVolumeIndexType maskIndex;
        
        // the seed is the index
        for( unsigned int col=0; col<3; col++ ) {
          maskIndex[ col ] = static_cast< int >( seedIndex[ col ] );
        }
        
        // if the pixel is masked, then we want to unmask it
        if( 0 == mask->GetPixel( maskIndex ) ) {
          mask->SetPixel( maskIndex, 1 );
          isUnmasked = true;
        }
      
      }
      
    }
  
  
  }
  
  if( isUnmasked ) {
    this->InvokeEvent( SeedsUnmaskedEvent() );  
  }

}


} // end namespace itk

#endif
