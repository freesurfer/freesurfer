#ifndef kvlTetrahedronInteriorConstIterator_h
#define kvlTetrahedronInteriorConstIterator_h

#include "itkImageConstIteratorWithIndex.h"
#include "kvlAtlasMesh.h"


namespace kvl
{
 
  
/**
 *
 * Iterator class that visits all the voxels of an image that lie inside a tetradron with
 * vertex coordinates p0, p1, p2, and p3, while also providing access to linearly interpolated
 * values at that location of sets of user-specified scalar values alpha0, alpha1, alpha2, alpha3 
 * at the vertices, i.e., 
 * 
 *   alpha = alpha0 * pi0 + alpha1 * pi1 + alpha2 * pi2 + alpha3 * pi3 
 * 
 * where pi0, pi1, pi2, and pi3 denote the baricentric coordinates of the voxel currently visited.
 * 
 * Typical usage therefore is something like this:
 * 
 *   TetrahedronInteriorIterator< ImageType::PixelType >  it( image, p0, p1, p2, p3 );
 *   it.AddExtraLoading( alpha0, alpha1, alpha2, alpha3 );
 *   for ( ; !it.IsAtEnd(); ++it )
 *     {
 *     std::cout << it.Value() << std::endl;
 *     std::cout << it.GetExtraLoadingInterpolatedValue( 0 ) << std::endl;
 *     }
 *
 * The baricentric coordinates are always available, even if the user doesn't ask for them --
 * this is because those baricentric coordinates are needed to decide whether or not a voxel
 * lies inside the tetradron. Internally this is accomplished by providing default "one-hot" 
 * loadings -- e.g., (alpha0,alpha1,alpha2,alpha3) = (0,1,0,0) will return second component of the 
 * baricentric coordinate. The baricentric coordinates can be accessed by it.GetPi0(), it.GetPi1(),
 * it.GetPi2(), and it.GetPi3().
 * 
 * At the core of the implementation is the fact that obtaining a linearly interpolated scalar value
 * "alpha" at a location "y" (3x1 vector of Eucledian coordinates) simply involves matrix multiplication:
 * 
 *   alpha( y ) = a * y + b
 * 
 * where a is 1x3 vector and b is a scalar.
 * This now means that if we know the interpolated value at one location, the interpolated value at a 
 * location "delta" (3x1 vector in Eucledian coordinates) away will be given by
 * 
 *   alpha( y + delta ) = a * (y+delta) + b = alpha( y ) + a * delta.
 * 
 * Now, because we know that voxels are set in a regular image grid, we now that visiting the voxel in 
 * the next row (i.e, delta = (1,0,0)) simply means adding a(1) (first component of the "a" row vector)
 * to whatever it is we currently have. Similarly, going to the next column adds a(2), and going to the
 * next slice a(3). 
 * 
 * Now, how do we calculate the values in a? (Note that we don't actually care for what's in b, as that's 
 * of no use to us.) First, we need to map the 3x1 vector of Eucledian coordinates 
 * "y" into baricentric ones (pi0,pi1,pi2,pi3). Given vertex coordinates p0, p1, p2, and p3, this
 * accomplished by doing
 *
 *   x = M * ( y - t );
 *   pi1 = x(1);
 *   pi2 = x(2);
 *   pi3 = x(3);
 *   pi0 = 1 - pi1 - pi2 - pi3;
 *
 * where 
 *
 *   M = inv( [ p1-p0 p2-p0 p3-p0 ] );
 *
 * and
 *
 *  t = p0;
 * 
 * To see why this is the case, consider the opposite direction: a tetradron with vertices
 * (0,0,0), (1,0,0), (0,1,0), and (0,0,1) will be mapped onto one with vertices p0, p1, p2, and p3
 * by doing
 * 
 * y =  [ p1-p0 p2-p0 p3-p0 ] * x + p0
 * 
 * Writing this in matrix form, we have that
 * 
 *  [ pi0 pi1 pi2 pi3 ]' = [ -1 -1 -1; 1 0 0; 0 1 0; 0 0 1 ] * x + [ 1 0 0 0 ]'
 *                       = [ -1 -1 -1; 1 0 0; 0 1 0; 0 0 1 ] * M * y + something
 * 
 * where I don't care about the 4x1 vector "something" as it will only affect "b" (in which I'm not
 * actually interested). This now means that 
 *
 *   alpha = [ alpha0 alpha1 alpha2 alpha3 ] * [ pi0 pi1 pi2 pi3 ]'
 *         = [ alpha0 alpha1 alpha2 alpha3 ] * [ -1 -1 -1; 1 0 0; 0 1 0; 0 0 1 ] * M * y + whatever
 *
 * so that finally 
 * 
 *   a = [ alpha0 alpha1 alpha2 alpha3 ] * [ -1 -1 -1; 1 0 0; 0 1 0; 0 0 1 ] * M.
 * 
 */  
  
  
// https://itk.org/Doxygen/html/classitk_1_1ImageConstIteratorWithIndex.html
  
template< typename TPixel >
class TetrahedronInteriorConstIterator : private itk::ImageConstIteratorWithIndex< typename itk::Image< TPixel, 3 > >
{
public:
#ifdef GEMS_DEBUG_RASTERIZE_VOXEL_COUNT
  int m_totalVoxel;
  int m_totalVoxelInTetrahedron;
#endif

  /** Standard class typedefs. */
  typedef TetrahedronInteriorConstIterator Self;
  typedef itk::ImageConstIteratorWithIndex< itk::Image< TPixel, 3 > > Superclass;

  /**
   * Index typedef support. While these were already typdef'ed in the superclass,
   * they need to be redone here for this subclass to compile properly with gcc.
   */
  /** Types inherited from the Superclass */
  typedef typename Superclass::IndexType             IndexType;
  typedef typename Superclass::SizeType              SizeType;
  typedef typename Superclass::OffsetType            OffsetType;
  typedef typename Superclass::RegionType            RegionType;
  typedef typename Superclass::ImageType             ImageType;
  typedef typename Superclass::PixelContainer        PixelContainer;
  typedef typename Superclass::PixelContainerPointer PixelContainerPointer;
  typedef typename Superclass::InternalPixelType     InternalPixelType;
  typedef typename Superclass::PixelType             PixelType;
  
  /** */
  typedef AtlasMesh::PointType   PointType;

  /** */
  typedef typename ImageType::OffsetValueType       OffsetValueType;

  /** Run-time type information (and related methods). */
  itkTypeMacroNoParent(TetrahedronInteriorConstIterator);

  /** Constructor */
  TetrahedronInteriorConstIterator( const ImageType *ptr,
                                    const PointType& p0, 
                                    const PointType& p1, 
                                    const PointType& p2, 
                                    const PointType& p3 );
  
  /** */
  const double& GetPi0() const
    {
    return m_InterpolatedValues[ 0 ];
    } 
  
  /** */
  const double& GetPi1() const
    {
    return m_InterpolatedValues[ 1 ];
    } 

  /** */
  const double& GetPi2() const
    {
    return m_InterpolatedValues[ 2 ];
    } 
    
  /** */
  const double& GetPi3() const
    {
    return m_InterpolatedValues[ 3 ];
    }
  
  /** Expose some of the basic iterator functionality we're going to use */
  using Superclass::IsAtEnd;
  using Superclass::GetIndex;
  using Superclass::Value;
  using Superclass::operator!=;
  
  /** Go to the next voxel that's inside the tetrahedron */                              
  Self&  operator++(); 
  
  /** */
  void AddExtraLoading( const double& alpha0, const double& alpha1, const double& alpha2, const double& alpha3 );
  
  /** */
  const double&  GetExtraLoadingInterpolatedValue( int extraLoadingNumber ) const
    {
    return m_InterpolatedValues[ 4 + extraLoadingNumber ];
    }
    
  /** */
  const double& GetExtraLoadingNextRowAddition( int extraLoadingNumber ) const
    {
    return m_NextRowAdditions[ 4 + extraLoadingNumber ];
    }

  /** */
  const double& GetExtraLoadingNextColumnAddition( int extraLoadingNumber ) const
    {
    return m_NextColumnAdditions[ 4 + extraLoadingNumber ];
    }
    
  /** */
  const double& GetExtraLoadingNextSliceAddition( int extraLoadingNumber ) const
    {
    return m_NextSliceAdditions[ 4 + extraLoadingNumber ];
    }

    
protected: //made protected so other iterators can access
  
  // 
#ifndef USING_STATIC_ARRAY
  std::vector< double >  m_InterpolatedValues;
  
  // 
  std::vector< double >  m_NextRowAdditions; 
  std::vector< double >  m_NextColumnAdditions; 
  std::vector< double >  m_NextSliceAdditions; 
#else
  //
  int m_NumberOfLoadings;
  
  double m_InterpolatedValues[MAX_LOADINGS];
  double m_NextRowAdditions[MAX_LOADINGS]; 
  double m_NextColumnAdditions[MAX_LOADINGS]; 
  double m_NextSliceAdditions[MAX_LOADINGS];
#endif
  
  // Make the data pointer visible to our subclasses
  using Superclass::m_Position;
  
private:

  //
  TetrahedronInteriorConstIterator(const Self &); // Not implemented
  void operator=(const Self &); // Not implemented
  
  // Go the next voxel inside the bounding box around the tetrahedron
  void MoveOnePixel(); 
  
  // Check if the current pixel is outside of the tetrahedron
  bool IsOutsideTetrahdron() const;
  
  //
  static bool CheckBorderCase( double a, double b, double c );

  // Things to help us backtrack to beginning of column and slice
#ifndef USING_STATIC_ARRAY
  std::vector< double >  m_ColumnBeginInterpolatedValues;
  std::vector< double >  m_SliceBeginInterpolatedValues;
#else
  double m_ColumnBeginInterpolatedValues[MAX_LOADINGS];
  double m_SliceBeginInterpolatedValues[MAX_LOADINGS];
#endif
  
  const InternalPixelType*  m_SliceBeginPosition;
  const InternalPixelType*  m_ColumnBeginPosition;

  
};

  
} // end namespace kvl

#include "kvlTetrahedronInteriorConstIterator.hxx"

#endif

