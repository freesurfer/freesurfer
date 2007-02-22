#ifndef __itkPoistatsFilter_h
#define __itkPoistatsFilter_h

#include <itkImageToImageFilter.h>
#include <itkArray2D.h>
#include <itkArray.h>

#include <itkScalarVector.h>

#include "dmri_poistats/datamodel/data/PoistatsReplicas.h"

#include "dmri_poistats/datamodel/PoistatsModel.h"

#include "dmri_poistats/datamodel/events/PoistatsEvents.h"

namespace itk
{

// the class should be named by algorithm-input-concept
//    poistats-image-filter

/** \class PoistatsFilter
 * Poistats is an itk filter for finding the most probable path between seed
 * regions.  Input a diffusion tensor image volume and seed regions, wait as the
 * replicas search the space and find the most probably path.
 */
template <class TInputImage, class TOutputImage>
class ITK_EXPORT PoistatsFilter :
  public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef PoistatsFilter                                     Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage>      Superclass;
  typedef SmartPointer<Self>                                 Pointer;
  typedef SmartPointer<const Self>                           ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(PoistatsFilter, ImageToImageFilter);
  
  /** Type for input image */
  typedef TInputImage InputImageType;

  typedef typename InputImageType::Pointer        InputImagePointer;
  typedef typename InputImageType::ConstPointer   InputImageConstPointer;

  /** Image typedef support. */
  typedef typename InputImageType::PixelType         PixelType;
  typedef typename InputImageType::RegionType        RegionType;
  typedef typename InputImageType::SizeType          SizeType;
  typedef typename InputImageType::IndexType         IndexType;
  typedef typename InputImageType::PointType         ContinuousIndexType;

  /** Type for output image */
  typedef TOutputImage OutputImageType;

  typedef typename OutputImageType::Pointer        OutputImagePointer;
  typedef typename OutputImageType::ConstPointer   OutputImageConstPointer;

  typedef typename OutputImageType::PixelType         OutputPixelType;
  typedef typename OutputImageType::RegionType        OutputRegionType;
  typedef typename OutputImageType::SizeType          OutputSizeType;
  typedef typename OutputImageType::IndexType         OutputIndexType;
  typedef typename OutputImageType::SpacingType       OutputSpacingType;
  
  /** Type for seed volume */
  typedef int SeedType;
  typedef itk::Image< SeedType, 3 > SeedVolumeType;
  typedef SeedVolumeType::Pointer SeedVolumePointer;
  typedef SeedVolumeType::PixelType SeedVolumePixelType;
  typedef SeedVolumeType::RegionType SeedVolumeRegionType;
  typedef SeedVolumeType::SizeType SeedVolumeSizeType;
  typedef SeedVolumeType::IndexType SeedVolumeIndexType;
    
  /** Type for sampling volume */
  typedef float SamplingType;
  typedef itk::Image< SamplingType, 3 > SamplingVolumeType;
  typedef SamplingVolumeType::Pointer SamplingVolumePointer;  
  typedef SamplingVolumeType::IndexType SamplingVolumeIndexType;  
  
  /** Type for mask volume */
  typedef int MaskType;
  typedef itk::Image< MaskType, 3 > MaskVolumeType;
  typedef MaskVolumeType::Pointer MaskVolumePointer;  
  typedef MaskVolumeType::PixelType MaskVolumePixelType;
  typedef MaskVolumeType::RegionType MaskVolumeRegionType;
  typedef MaskVolumeType::SizeType MaskVolumeSizeType;
  typedef MaskVolumeType::IndexType MaskVolumeIndexType;
    
  typedef itk::Array2D< double > MatrixType;
  typedef MatrixType* MatrixPointer;
  typedef std::vector< MatrixPointer > MatrixListType;
  
  typedef itk::Array< double > ArrayType;
  typedef ArrayType* ArrayPointer;
  typedef std::vector< ArrayPointer > ArrayListType;
  
  typedef itk::Array< int > VoxelIndexType;
  typedef VoxelIndexType* VoxelIndexPointer;
  
  /** Type for ODF look up table */
  typedef itk::Image< int, 3 > OdfLookUpTableType;  
  typedef OdfLookUpTableType::Pointer OdfLookUpTablePointer;
  typedef OdfLookUpTableType::SizeType OdfLookUpSizeType;
  typedef OdfLookUpTableType::IndexType OdfLookUpIndexType;
  typedef OdfLookUpTableType::RegionType OdfLookUpRegionType;

  // used with get and set output
  itkStaticConstMacro( PATH_DENSITY_OUTPUT, int, 0 );
  itkStaticConstMacro( OPTIMAL_PATH_DENSITY_OUTPUT, int, 1 );

  itkStaticConstMacro( INVALID_INDEX, int, -1 );
  
  itkGetMacro( InitialSigma, int );
  itkSetMacro( InitialSigma, int );

  itkGetMacro( NumberOfSamplePoints, int );
  itkSetMacro( NumberOfSamplePoints, int );

  OutputImageType *GetOptimalDensity();
  
  OutputImageType *GetDensity();
  
  void CalculateTensor2Odf( itk::Matrix< double, 3, 3 > *tensor,
    ArrayPointer odf );
    
  double CalculateOdfPathEnergy( itk::Array2D< double > *path,
    itk::Array< double > **odfs, ArrayPointer outputEnergies );
        
  static void RoundPath( itk::Array2D< int > *outputRoundedArray, 
                         itk::Array2D< double >* inputArray );
        
  itkGetConstMacro( Exchanges, int );

  itkGetConstMacro( CurrentIteration, int );
  itkGetConstMacro( CurrentLull, int );
  itkGetConstMacro( CurrentEnergyDifference, double );
  itkGetConstMacro( GlobalMinEnergy, double );
  itkGetConstMacro( ElapsedTime, double );  
  
  MatrixPointer GetStartSeeds();
  MatrixPointer GetEndSeeds();

  void SetRandomSeed( const long seed );
  
  void GetOdfsAtPoints( itk::Array< double >** outputOdfs, 
                        itk::Array2D< int >* inputPoints );

  // TODO: this method creates new odfs without deallocating the memory later
  void ConstructOdfList();
    
  static double GetDistance( vnl_vector< double >* point1, 
                             itk::Array< double >* point2 );
                        
  static void GetPointClosestToCenter( itk::Array2D< double >* seeds, 
                                       itk::Array< double >* closestSeed );

  static void GetCenterOfMass( itk::Array2D< double >* mass, 
                               itk::Array< double >* center );
    
  void GetInitialPoints( itk::Array2D< double >* initialPoints );
  
  int GetNumberOfInitialPoints() const;  
  
  itkGetMacro( Polarity, MatrixType );
  itkSetMacro( Polarity, MatrixType );
  
  itkGetMacro( FinalPath, MatrixType );
  itkSetMacro( FinalPath, MatrixType );
  
  void ConvertPointsToImage( MatrixPointer points, OutputImagePointer volume );
  
  static void GetPositiveMinimumInt( MatrixPointer points, const int radius, 
    VoxelIndexPointer minimum );

  static void GetPositiveMaximumInt( MatrixPointer points, const int radius, 
    VoxelIndexPointer maximum, OutputRegionType rowCeiling );
    
  void SetSeedVolume( SeedVolumePointer volume );

  void GetPathProbabilities( MatrixPointer path, ArrayPointer probabilities );

  int GetNumberOfSteps();
  void SetNumberOfSteps( const int nSteps );

  itkGetMacro(NumberOfDirections, int);
  itkSetMacro(NumberOfDirections, int);

  void GetAggregateReplicaDensities(
    MatrixListType replicaPaths, 
    OutputImagePointer aggregateDensities );

  itkGetMacro(MaxTime, int);
  itkSetMacro(MaxTime, int);
  
  void CalculateFinalPathProbabilities();  

  void CalculateBestPathProbabilities();

  MatrixType GetBestPathsProbabilities();
  
  int GetNumberOfReplicas();
  void SetNumberOfReplicas( const int nReplicas );

  itkGetMacro(FinalPathProbabilities, ArrayType);
//  itkSetMacro(FinalPathProbabilities, ArrayType);

  itkSetMacro( SamplingVolume, SamplingVolumePointer );
//  void SetSamplingVolume( SamplingVolumePointer volume );
  ArrayType GetSamples();
  
  itkGetMacro( MaskVolume, MaskVolumePointer );
  void SetMaskVolume( MaskVolumePointer volume );
  
  OdfLookUpTablePointer GetOdfLookUpTable();  
  
  void SetNextSeedValueToUse( int seedValue );
  
  static void GenerateRotationMatrix3u2v( 
    ArrayPointer u, ArrayPointer v, MatrixPointer outputRotationMatrix );
    
  void GetMagnetToSliceFrameRotation( MatrixPointer rotation );

  void ParseSeedVolume();

  static void CalculateDensityMatrix( itk::Array2D< double > *angles, 
    itk::Array2D< double > *densityMatrix );

  itkGetMacro(MaxLull, int);
  itkSetMacro(MaxLull, int);
  
  int GetNumberOfControlPoints();
  void SetNumberOfControlPoints( const int nPoints );  
  
  itkGetMacro( ReplicaExchangeProbability, double );
  itkSetMacro( ReplicaExchangeProbability, double );  
  
  itkGetMacro( SigmaTimeConstant, double );
  itkSetMacro( SigmaTimeConstant, double );

  itkGetMacro( PointsToImageGamma, double );
  itkSetMacro( PointsToImageGamma, double );
      
  double GetCurrentMeanOfEnergies() const;

  double GetCurrentMinOfEnergies() const;      

protected:
  PoistatsFilter();
  virtual ~PoistatsFilter();
  
  virtual void GenerateOutputInformation(){}; // do nothing
  virtual void GenerateData();  

private:
  static const double NO_ZERO_SHELL_252[][3];
  
  PoistatsModel *m_PoistatsModel;

  SeedVolumePointer m_SeedVolume;
  
  SamplingVolumePointer m_SamplingVolume;

  MaskVolumePointer m_MaskVolume;
  
  OdfLookUpTablePointer m_OdfLookUpTable;
   
  ArrayListType m_Odfs;

  int m_InitialSigma;
  
  itkStaticConstMacro( DEFAULT_NUMBER_OF_SAMPLE_POINTS, int, 100 );

  int m_NumberOfSamplePoints;
  
  itkStaticConstMacro( DEFAULT_NUMBER_OF_REPLICAS, int, 100 );
  itkStaticConstMacro( DEFAULT_NUMBER_OF_STEPS, int, 50 );
  
  itkStaticConstMacro( DEFAULT_MAX_LULL, int, 10 );
  itkStaticConstMacro( DEFAULT_MAX_TIME, int, 300 );
  int m_MaxLull; 
  int m_MaxTime; 
  
  int m_Exchanges;
  
  itkStaticConstMacro( DEFAULT_NUMBER_OF_DIRECTIONS, int, 252 );
  int m_NumberOfDirections;

  /* itkStaticConstMacro can only be used with ints because:
   The StaticConstMacro was designed for dealing with numeric
   parameters of templated classes. Such parameter can only
   be "int", "unsigned int", and "enum".
  
   However, some compilers do not support the enums, while other
   do not support the initialization of static variables.
  
   You could be ok with using the StaticConstMacro with a double
   type *ONLY* in some compilers... 
 */
  // itkStaticConstMacro( DEFAULT_COOL_FACTOR, double, 0.995 );
  static const double DEFAULT_COOL_FACTOR = 0.995;

  double m_CoolFactor;
  itkGetMacro(CoolFactor, double);
  itkSetMacro(CoolFactor, double);
  
  std::vector< MatrixPointer > m_Seeds;

  std::vector< int > m_SeedValuesToUse;
  
  MatrixType m_Polarity;
  
  MatrixType m_FinalPath;
  
  void InitPaths();
  
  itkSetMacro(Exchanges, int);
  
  ArrayType m_SliceUp;
                                                        
  static void CalculateAnglesBetweenVectors( itk::Array2D< double > *vectors, 
    itk::Array< double > *angles );

  bool ShouldUpdateEnergy( double currentEnergy, double previousEnergy, 
                           double currentTemperature );
                           
  void AllocateOutputImage(  const int outputIndex, OutputImagePointer image);
  
  static void GetSortedUniqueSeedValues( SeedVolumePointer volume, 
    std::vector< std::pair< SeedType, int > > *seedValues );
    
  ArrayType m_FinalPathProbabilities;

  static void TakeUnionOfSeeds(
    std::vector< std::pair< SeedType, int > > *seeds1,
    std::vector< SeedType > *seeds2 );
  
  typedef vnl_matrix< double > VnlMatrixType;
    
  VnlMatrixType m_TensorGeometry;  
  VnlMatrixType GetTensorGeometry();
  
  typedef vnl_vector< double > VnlVectorType;
  ArrayType m_InvalidOdf;
      
  PoistatsReplicas *m_Replicas;
  
  double m_ReplicaExchangeProbability;
  static const double DEFAULT_REPLICA_EXCHANGE_PROBABILITY = 0.05;

  double m_SigmaTimeConstant;
  static const double DEFAULT_SIGMA_TIME_CONSTANT = 2*1e2;  

  double m_PointsToImageGamma;
  static const double DEFAULT_POINTS_TO_IMAGE_GAMMA = 0.5;  
  
  int m_CurrentIteration;
  int m_CurrentLull;

  double m_GlobalMinEnergy;
  itkSetMacro( GlobalMinEnergy, double );

  double m_ElapsedTime;
  itkSetMacro( ElapsedTime, double );  

  double m_CurrentEnergyDifference;
  
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPoistatsFilter.txx"
#endif

#endif
