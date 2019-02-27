#ifndef __PoistatsModel_h
#define __PoistatsModel_h

#include "bspline/itkBSplineDataPointSetToImageFilter.h" // itkutils

#include <vnl/vnl_random.h>
#include <vector>

#include <itkArray2D.h>
#include <itkArray.h>
#include <itkPointSet.h>

 
#include "mri.h"


class InitializePath;

class PoistatsModel
{
public:

  enum {
    INITIALIZE_NORMAL = 0,
    INITIALIZE_EIGEN_VECTOR = 1,
    INITIALIZE_FIELD_LINE = 2,
  };

  typedef itk::Array2D< double > MatrixType;
  typedef MatrixType* MatrixPointer;

  typedef itk::Array< double > ArrayType;
  typedef ArrayType* ArrayPointer;

  // cubic spline stuff
  typedef itk::Vector< double, 3 >    VectorType;
  typedef itk::Image< VectorType, 1 > OutputImageType;
  typedef itk::PointSet< VectorType, 1 > PointSetType;
  typedef itk::BSplineDataPointSetToImageFilter
     < PointSetType, OutputImageType >  CubicSplineFilterType;
  typedef CubicSplineFilterType::Pointer CubicSplineFilterPointer;

  PoistatsModel();
  PoistatsModel( const long seed );
  ~PoistatsModel();

  void SetRandomSeed( const long seed );
      
  double GetRandomNumber();
  
  double GetNormallyDistributedRandomNumber();

  int GetRandomInt( const int floor, const int ceiling );
  
  double GetRandomNumberWithoutRange();
  
  bool IsDebug();

  void SetDebug( const bool isDebug );
  
  bool IsUsingPathInitialization();

  /**
   * Returns the eigenvectors.
   */
  MRI* GetEigenVectors();
  void SetEigenVectors( MRI* vectors );   
  
  /**
   * Returns the volume with the seed labels.
   */
  MRI* GetSeedVolume();
  void SetSeedVolume( MRI* volume );   
  
  /**
   * Returns the seed values to be used.
   */
  std::vector< int >* GetSeedValues();
  void SetSeedValues( std::vector< int >* values);

  MatrixPointer RethreadPath( 
    MatrixPointer originalPath, const int nNewSamples );

  MatrixPointer CubicSplineInterpolation( 
    MatrixPointer originalPath,
    ArrayPointer originalPathGrid,
    const int nNewSamples );

  void SetNumberOfControlPoints( const int nPoints );

  int GetNumberOfControlPoints();

  static void SpaceEvenly( ArrayPointer outputArray, const double floor, 
    const double ceiling );
    
  static void CalculatePathMagnitude( 
    MatrixPointer path, ArrayPointer magnitude );

  static void CalculateCumulativeSum(
    ArrayPointer inputArray, ArrayPointer cumulativeSum );    
    
  static void CalculatePathVectors( MatrixPointer path, MatrixPointer vectors );

  static void CalculateMagnitude( MatrixPointer path, ArrayPointer magnitude );
  
  int GetInitialSigma() const;

  void SetInitialSigma( const int initialSigma );
  
  InitializePath *GetPathInitializer() const;
  
  void SetInitializePathMode( const int mode );

  int GetInitializePathMode() const;
  
  void SetFieldLineRadius( const double radius );
  
  double GetFieldLineRadius() const;
      
private:

  vnl_random m_RandomNumberGenerator;
  
  bool m_IsDebug;
  
  MRI* m_EigenVectors;
  MRI* m_SeedVolume;
  std::vector< int >* m_SeedValues;

  CubicSplineFilterPointer m_CubicSplineFilter;

  int m_NumberOfControlPoints;
  
  int m_InitialSigma;
  
  InitializePath *m_PathInitializer;
  
  int m_InitializeMode;
  
  double m_FieldLineRadius;
    
  /**
   * Initializes the class.
   */
  void Init();
  
  void FreeVector( std::vector< int > *v );  
  
};

#endif

