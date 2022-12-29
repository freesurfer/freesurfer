#ifndef kvlTetrahedronInteriorConstIterator_hxx
#define kvlTetrahedronInteriorConstIterator_hxx

#include "kvlTetrahedronInteriorConstIterator.h"

namespace kvl
{


//
//
//  
template< typename TPixel >
TetrahedronInteriorConstIterator< TPixel >
::TetrahedronInteriorConstIterator( const ImageType *ptr,
                                    const PointType& p0, 
                                    const PointType& p1, 
                                    const PointType& p2, 
                                    const PointType& p3 )
#ifndef USING_STATIC_ARRAY
: Superclass( ptr, RegionType() ), 
  m_InterpolatedValues( 4 ), 
  m_NextRowAdditions( 4 ), 
  m_NextColumnAdditions( 4 ), 
  m_NextSliceAdditions( 4 ),
  m_ColumnBeginInterpolatedValues( 4 ),
  m_SliceBeginInterpolatedValues( 4 )
#else
: Superclass()  //Superclass( ptr, RegionType() )
#endif
{
  
#ifdef GEMS_DEBUG_RASTERIZE_VOXEL_COUNT
  m_totalVoxel = 0;
  m_totalVoxelInTetrahedron = 0;
#endif

  // ============================================================================================
  //
  // Part I: Compute a valid bounding box around the tethradron specified by vertices (p0,p1,p2,p3). 
  // If the tetradron falls outside the buffered image region, the bounding box is clipped accordingly
  //
  // ============================================================================================
  typedef typename ImageType::RegionType   RegionType;
  typedef typename RegionType::IndexType  IndexType;
  typedef typename IndexType::IndexValueType  IndexValueType;
  typedef typename RegionType::SizeType  SizeType;
  typedef typename SizeType::SizeValueType  SizeValueType;
  
  // Compute the coordinates of the lower corner of the bounding box around the tetradron
  PointType  lowerCorner = p0;
  for ( int i = 0; i < 3; i++ )
    {
    if ( p1[ i ] < lowerCorner[ i ] )
      {
      lowerCorner[ i ] = p1[ i ];  
      }  
    if ( p2[ i ] < lowerCorner[ i ] )
      {
      lowerCorner[ i ] = p2[ i ];  
      }  
    if ( p3[ i ] < lowerCorner[ i ] )
      {
      lowerCorner[ i ] = p3[ i ];  
      }  
    }
  //std::cout << "lowerCorner: " << lowerCorner << std::endl;  

  // Compute the coordinates of the upper corner of the bounding box around the tetradron
  PointType  upperCorner = p0;
  for ( int i = 0; i < 3; i++ )
    {
    if ( p1[ i ] > upperCorner[ i ] )
      {
      upperCorner[ i ] = p1[ i ];  
      }  
    if ( p2[ i ] > upperCorner[ i ] )
      {
      upperCorner[ i ] = p2[ i ];  
      }  
    if ( p3[ i ] > upperCorner[ i ] )
      {
      upperCorner[ i ] = p3[ i ];  
      }
    
    }
  //std::cout << "upperCorner: " << upperCorner << std::endl;  
    
  // Compute the lower cornder index, while clipping to the buffered region
  IndexType  lowerCornerIndex = ptr->GetBufferedRegion().GetIndex();
  for ( int i = 0; i < 3; i++ )
    {
    if ( lowerCorner[ i ] > lowerCornerIndex[ i ] )
      {
      lowerCornerIndex[ i ] = itk::Math::Ceil< IndexValueType >( lowerCorner[ i ] );
      }
      
    // Pathological case where tethradron is completely outside of image domain;
    // let's make sure the size of our region is then 0
    if ( lowerCornerIndex[ i ] > ptr->GetBufferedRegion().GetUpperIndex()[ i ] )
      {
      lowerCornerIndex[ i ] = ptr->GetBufferedRegion().GetUpperIndex()[ i ] + 1;
      }  
    }
  //std::cout << "lowerCornerIndex: " << lowerCornerIndex << std::endl;
  

  // Compute the upper cornder index, while clipping to the buffered region
  IndexType  upperCornerIndex = ptr->GetBufferedRegion().GetUpperIndex();
  for ( int i = 0; i < 3; i++ )
    {
    if ( upperCorner[ i ] < upperCornerIndex[ i ] )
      {
      upperCornerIndex[ i ] = itk::Math::Floor< IndexValueType >( upperCorner[ i ] );
      }
      
    // Pathological case where tethradron is completely outside of image domain;
    // let's make sure the size of our region is then 0
    if ( upperCornerIndex[ i ] < ptr->GetBufferedRegion().GetIndex()[ i ] )
      {
      upperCornerIndex[ i ] = ptr->GetBufferedRegion().GetIndex()[ i ] - 1;
      }  
    }
  //std::cout << "upperCornerIndex: " << upperCornerIndex << std::endl;
  

  // ============================================================================================
  //
  // Part II: Set up the base class image iterator stuff.
  //
  // ============================================================================================
  RegionType region;
  region.SetIndex( lowerCornerIndex );
  region.SetUpperIndex( upperCornerIndex );
  //std::cout << "region: " << region << std::endl;
  Superclass::operator=( Superclass( ptr, region ) ); // Workaround for non-existing this->SetRegion( region )
  
  m_SliceBeginPosition = this->m_Position;
  m_ColumnBeginPosition = this->m_Position;
  
  
  // ============================================================================================
  //
  // Part III: Precompute some matrices that will allow us to map a 3x1 vector of Eucledian coordinates 
  // ("y") into baricentric ones (pi0,pi1,pi2,pi3). Given vertex coordinates p0, p1, p2, and p3, this
  // accomplished by doing
  //
  //   x = M * ( y - t );
  //   pi1 = x(1);
  //   pi2 = x(2);
  //   pi3 = x(3);
  //   pi0 = 1 - pi1 - pi2 - pi3;
  //
  // where 
  //
  //   M = inv( [ p1-p0 p2-p0 p3-p0 ] );
  //
  // and
  //
  //   t = p0;
  //
  // To see why this is the case, consider the opposite direction: a tetradron with vertices
  // (0,0,0), (1,0,0), (0,1,0), and (0,0,1) will be mapped onto one with vertices p0, p1, p2, and p3
  // by doing
  //
  //  y =  [ p1-p0 p2-p0 p3-p0 ] * x + p0;  
  //
  // ============================================================================================
  
  // t = p0
  const double  t1 = p0[ 0 ];
  const double  t2 = p0[ 1 ];
  const double  t3 = p0[ 2 ];
  
  // M = inv( [ p1-p0 p2-p0 p3-p0 ] )
  // where the inversion of a 3x3 matrix is given by:
  //
  // (cf.  https://en.wikipedia.org/wiki/Invertible_matrix#Inversion_of_3.C3.973_matrices)
  const double  a = p1[0] - p0[0];
  const double  b = p2[0] - p0[0];
  const double  c = p3[0] - p0[0];
  const double  d = p1[1] - p0[1];
  const double  e = p2[1] - p0[1];
  const double  f = p3[1] - p0[1];
  const double  g = p1[2] - p0[2];
  const double  h = p2[2] - p0[2];
  const double  i = p3[2] - p0[2];
  
  const double  A = ( e * i - f * h );
  const double  D = -( b * i - c * h );
  const double  G = ( b * f - c * e );
  const double  B = -(d * i - f * g );
  const double  E = ( a * i - c * g );
  const double  H = -( a * f - c * d );
  const double  C = ( d * h - e * g );
  const double  F = - (a * h - b * g );
  const double  I = ( a * e - b * d );
  
  const double  determinant = a * A + b * B + c * C;
  // M = 1/determinant * [ A D G; B E H; C F I ]
  const double  m11 = A / determinant;
  const double  m21 = B / determinant;
  const double  m31 = C / determinant;
  const double  m12 = D / determinant;
  const double  m22 = E / determinant;
  const double  m32 = F / determinant;
  const double  m13 = G / determinant;
  const double  m23 = H / determinant;
  const double  m33 = I / determinant;
  
  
  // ============================================================================================
  //
  // Part IV: Precompute the baricentric coordinates of the first voxel of the first column in 
  // the first slice (i.e, the one we're currently pointing to)
  //
  // ============================================================================================
  const double  YminT1 = this->GetIndex()[ 0 ] - t1;
  const double  YminT2 = this->GetIndex()[ 1 ] - t2;
  const double  YminT3 = this->GetIndex()[ 2 ] - t3;
  
  // 
  const double  pi1 = m11 * YminT1 + m12 * YminT2 + m13 * YminT3;
  const double  pi2 = m21 * YminT1 + m22 * YminT2 + m23 * YminT3;
  const double  pi3 = m31 * YminT1 + m32 * YminT2 + m33 * YminT3;
  const double  pi0 = 1.0 - pi1 - pi2 - pi3;

  //
  m_InterpolatedValues[ 0 ] = pi0;
  m_InterpolatedValues[ 1 ] = pi1;
  m_InterpolatedValues[ 2 ] = pi2;
  m_InterpolatedValues[ 3 ] = pi3;
 
#ifdef USING_STATIC_ARRAY
  m_NumberOfLoadings = 4;
#endif

  //
  for ( int loadingNumber = 0; loadingNumber < 4; loadingNumber++ ) 
    {
    m_ColumnBeginInterpolatedValues[ loadingNumber ] = m_InterpolatedValues[ loadingNumber ];
    m_SliceBeginInterpolatedValues[ loadingNumber ] = m_InterpolatedValues[ loadingNumber ];
    }
    
  //
  // m_NextRowAdditions, m_NextColumnAdditions, m_NextSliceAdditions are constants.
  // they are used to increment m_InterpolatedValues
  // [ -1 -1 -1; 1 0 0; 0 1 0; 0 0 1 ] * M * [1 0 0]' 
  m_NextRowAdditions[ 0 ] = -( m11 + m21 + m31 );
  m_NextRowAdditions[ 1 ] = m11;
  m_NextRowAdditions[ 2 ] = m21;
  m_NextRowAdditions[ 3 ] = m31;

  //
  // [ -1 -1 -1; 1 0 0; 0 1 0; 0 0 1 ] * M * [0 1 0]'
  m_NextColumnAdditions[ 0 ] = -( m12 + m22 + m32 );
  m_NextColumnAdditions[ 1 ] = m12;
  m_NextColumnAdditions[ 2 ] = m22;
  m_NextColumnAdditions[ 3 ] = m32;

  //
  // [ -1 -1 -1; 1 0 0; 0 1 0; 0 0 1 ] * M * [0 0 1]'
  m_NextSliceAdditions[ 0 ] = -( m13 + m23 + m33 );
  m_NextSliceAdditions[ 1 ] = m13;
  m_NextSliceAdditions[ 2 ] = m23;
  m_NextSliceAdditions[ 3 ] = m33;
  
  
  // ============================================================================================
  //
  // Part V: Advance to the first voxel that is actually inside the tetradron
  //
  // ============================================================================================
  while ( this->IsOutsideTetrahdron() && !this->IsAtEnd() )
    {
    this->MoveOnePixel();
    }
 
  
}


//
//
//
template< typename TPixel >
void
TetrahedronInteriorConstIterator< TPixel >
::AddExtraLoading( const double& alpha0, const double& alpha1, const double& alpha2, const double& alpha3 )
{
  
#ifndef USING_STATIC_ARRAY
  //
  m_InterpolatedValues.push_back( alpha0 * m_InterpolatedValues[ 0 ] + 
                                  alpha1 * m_InterpolatedValues[ 1 ] + 
                                  alpha2 * m_InterpolatedValues[ 2 ] +
                                  alpha3 * m_InterpolatedValues[ 3 ] );
  m_ColumnBeginInterpolatedValues.push_back( alpha0 * m_ColumnBeginInterpolatedValues[ 0 ] + 
                                             alpha1 * m_ColumnBeginInterpolatedValues[ 1 ] + 
                                             alpha2 * m_ColumnBeginInterpolatedValues[ 2 ] +
                                             alpha3 * m_ColumnBeginInterpolatedValues[ 3 ] );
  m_SliceBeginInterpolatedValues.push_back( alpha0 * m_SliceBeginInterpolatedValues[ 0 ] + 
                                            alpha1 * m_SliceBeginInterpolatedValues[ 1 ] + 
                                            alpha2 * m_SliceBeginInterpolatedValues[ 2 ] +
                                            alpha3 * m_SliceBeginInterpolatedValues[ 3 ] );

  m_NextRowAdditions.push_back( alpha0 * m_NextRowAdditions[ 0 ] + 
                                alpha1 * m_NextRowAdditions[ 1 ] + 
                                alpha2 * m_NextRowAdditions[ 2 ] +
                                alpha3 * m_NextRowAdditions[ 3 ] ); 
  m_NextColumnAdditions.push_back( alpha0 * m_NextColumnAdditions[ 0 ] + 
                                   alpha1 * m_NextColumnAdditions[ 1 ] + 
                                   alpha2 * m_NextColumnAdditions[ 2 ] +
                                   alpha3 * m_NextColumnAdditions[ 3 ] ); 
  m_NextSliceAdditions.push_back( alpha0 * m_NextSliceAdditions[ 0 ] + 
                                  alpha1 * m_NextSliceAdditions[ 1 ] + 
                                  alpha2 * m_NextSliceAdditions[ 2 ] +
                                  alpha3 * m_NextSliceAdditions[ 3 ] );
#else
  //
  m_InterpolatedValues[m_NumberOfLoadings] = alpha0 * m_InterpolatedValues[ 0 ] + 
                                             alpha1 * m_InterpolatedValues[ 1 ] + 
                                             alpha2 * m_InterpolatedValues[ 2 ] +
                                             alpha3 * m_InterpolatedValues[ 3 ];
  m_ColumnBeginInterpolatedValues[m_NumberOfLoadings] = alpha0 * m_ColumnBeginInterpolatedValues[ 0 ] + 
                                                        alpha1 * m_ColumnBeginInterpolatedValues[ 1 ] + 
                                                        alpha2 * m_ColumnBeginInterpolatedValues[ 2 ] +
                                                        alpha3 * m_ColumnBeginInterpolatedValues[ 3 ];
  m_SliceBeginInterpolatedValues[m_NumberOfLoadings] = alpha0 * m_SliceBeginInterpolatedValues[ 0 ] + 
                                                       alpha1 * m_SliceBeginInterpolatedValues[ 1 ] + 
                                                       alpha2 * m_SliceBeginInterpolatedValues[ 2 ] +
                                                       alpha3 * m_SliceBeginInterpolatedValues[ 3 ];

  m_NextRowAdditions[m_NumberOfLoadings] = alpha0 * m_NextRowAdditions[ 0 ] + 
                                           alpha1 * m_NextRowAdditions[ 1 ] + 
                                           alpha2 * m_NextRowAdditions[ 2 ] +
                                           alpha3 * m_NextRowAdditions[ 3 ]; 
  m_NextColumnAdditions[m_NumberOfLoadings] = alpha0 * m_NextColumnAdditions[ 0 ] + 
                                              alpha1 * m_NextColumnAdditions[ 1 ] + 
                                              alpha2 * m_NextColumnAdditions[ 2 ] +
                                              alpha3 * m_NextColumnAdditions[ 3 ]; 
  m_NextSliceAdditions[m_NumberOfLoadings] = alpha0 * m_NextSliceAdditions[ 0 ] + 
                                             alpha1 * m_NextSliceAdditions[ 1 ] + 
                                             alpha2 * m_NextSliceAdditions[ 2 ] +
                                             alpha3 * m_NextSliceAdditions[ 3 ];
  m_NumberOfLoadings++;
#endif
  
 
}  

//
//
//
template< typename TPixel >
TetrahedronInteriorConstIterator< TPixel >&  
TetrahedronInteriorConstIterator< TPixel >
::operator++()
{

  this->MoveOnePixel();
  while ( this->IsOutsideTetrahdron() && !this->IsAtEnd() )
    {
    this->MoveOnePixel();
    }

  return *this;
}

   
   
   
//
//
//
template< typename TPixel >
void
TetrahedronInteriorConstIterator< TPixel >
::MoveOnePixel()
{
  
#ifdef GEMS_DEBUG_RASTERIZE_VOXEL_COUNT
  m_totalVoxel++;
#endif

  // In principle, all we need to do is to do 
  //
  //   y = this->GetIndex()
  //
  // and then compute the baricentric coordinates using
  // 
  //   x = M * ( y - t )
  //
  // However, since our y's are nicely ordered on an image grid, we can save a lot of
  // multiplications by noticing that
  // 
  //  M * ( ( y + delta ) - t ) = M * ( y - t ) + M * delta
  //
  // where delta is typically (1,0,0) (so that M * delta is just the first column of M),
  // although sometimes (when we're at the end of our span) it will be (0,1,0) (second column of M),
  // and even more less frequently (0,0,1) (third column of M)

#ifndef USING_STATIC_ARRAY
  const int  numberOfLoadings = m_InterpolatedValues.size();
#else
  const int  numberOfLoadings = m_NumberOfLoadings;
#endif
  if ( this->m_PositionIndex[ 0 ] < ( this->m_EndIndex[ 0 ] - 1 ) )
    {
    // We can simply walk to the next row
    
    // Update the index
    this->m_PositionIndex[ 0 ]++;  
      
    // Update the data pointer
    this->m_Position++;  
 
    //  Update the baricentric coordinates
    // for ( int loadingNumber = 0; loadingNumber < numberOfLoadings; loadingNumber++ )
    //   {
    //   m_InterpolatedValues[ loadingNumber ] += m_NextRowAdditions[ loadingNumber ];
    //   }
    for ( int loadingNumber = 0; loadingNumber < numberOfLoadings; loadingNumber++ )
      {
      m_InterpolatedValues[ loadingNumber ] += m_NextRowAdditions[ loadingNumber ];  
      }  
  
    }  
  else if ( this->m_PositionIndex[ 1 ] < ( this->m_EndIndex[ 1 ] - 1 ) )
    {
    // We've reached the last row of our column. and we can simply go to 
    // the next column
  
    // Update the index
    this->m_PositionIndex[ 0 ] = this->m_BeginIndex[ 0 ];
    this->m_PositionIndex[ 1 ]++;  
      
    // Update the data pointer
    m_ColumnBeginPosition += this->m_OffsetTable[ 1 ];
    this->m_Position  =  m_ColumnBeginPosition;  
 
    //  Update the baricentric coordinates
    for ( int loadingNumber = 0; loadingNumber < numberOfLoadings; loadingNumber++ )
      {
      // m_ColumnBeginInterpolatedValues is the new base 
      m_ColumnBeginInterpolatedValues[ loadingNumber ] += m_NextColumnAdditions[ loadingNumber ];  
      m_InterpolatedValues[ loadingNumber ] = m_ColumnBeginInterpolatedValues[ loadingNumber ];
      }  

    
    }
  else if ( this->m_PositionIndex[ 2 ] < ( this->m_EndIndex[ 2 ] - 1 ) )
    {
    // We've reached the last row of the last column of our slice, but 
    // we can simply go to the next slice
      
    // Update the index
    this->m_PositionIndex[ 0 ] = this->m_BeginIndex[ 0 ];
    this->m_PositionIndex[ 1 ] = this->m_BeginIndex[ 1 ];
    this->m_PositionIndex[ 2 ]++;  
      
    // Update the data pointer
    m_SliceBeginPosition += this->m_OffsetTable[ 2 ];
    m_ColumnBeginPosition = m_SliceBeginPosition;  
    this->m_Position = m_SliceBeginPosition;  
 
    //  Update the baricentric coordinates
    for ( int loadingNumber = 0; loadingNumber < numberOfLoadings; loadingNumber++ )
      {
      // m_SliceBeginInterpolatedValues is the new base 
      // update m_ColumnBeginInterpolatedValues for next column calculation
      m_SliceBeginInterpolatedValues[ loadingNumber ] += m_NextSliceAdditions[ loadingNumber ];
      m_ColumnBeginInterpolatedValues[ loadingNumber ] = m_SliceBeginInterpolatedValues[ loadingNumber ];  
      m_InterpolatedValues[ loadingNumber ] = m_SliceBeginInterpolatedValues[ loadingNumber ];
      }
      
    }
  else  
    {
    this->m_Remaining = false;
    }  
    
}




//
//
//
template< typename TPixel >
bool
TetrahedronInteriorConstIterator< TPixel >
::IsOutsideTetrahdron() const
{
  
  // In general, a pixel falls outside the tetradron if one of its baricentric
  // coordinates turns negative. However this still leaves the hairy issue of
  // what to do with those special snowflake cases where one or more of the 
  // coordinates is *exactly* zero (meaning it lies on the face of one of the
  // four triangles composing the tetradron): these cases should only belong
  // to a single tetradron, as otherwise they will be visited by all tetradra
  // sharing the same face, which means they'll be counted multiple times when
  // evaluating e.g., a cost function that is the sum over all voxels
  //
  // Similar to the "top-left rule" for rasterizing triangles, we can come up
  // with similar rules for tetrahedra. The philosophy is that if a baricentric
  // coordinate is exactly zero, we (virtually) shift the voxel a tiny fraction
  // to the right; if this changes the baricentric coordinate to become positive,
  // the voxel will be considered inside. This will happen if the corresponding
  // element in the first column of the M matrix is positive. Of course it is 
  // possible that this elememt is exactly zero; if this is the case (meaning 
  // that we're lying on a face that is parellell with the x-axis), we can test
  // for the second direction (second column of M) -- virtually pushing the point,
  // "up"; and if also that element is zero (which means the face is perpendicular
  // to the z-axis) we look at the element in the third column of M (virtually
  // pushing the point along the z-axis by a tiny fraction)
  //
  // For a nice drawing of the "top-left rule" for triangles, see
  //   https://msdn.microsoft.com/en-us/library/windows/desktop/cc627092%28v=vs.85%29.aspx#Triangle
  // For a thorough exaplanation of the art of rasterizing triangles, see
  //   https://fgiesen.wordpress.com/2013/02/08/triangle-rasterization-in-practice/
  
  //const double&  pi0 = this->GetPi0();
  //const double&  pi1 = this->GetPi1();
  //const double&  pi2 = this->GetPi2();
  //const double&  pi3 = this->GetPi3();
  
  // Obvious culling first
  //if ( ( pi0 < 0 ) || ( pi1 < 0 ) || ( pi2 < 0 ) || ( pi3 < 0 ) )
  if ( ( m_InterpolatedValues[ 0 ] < 0 ) || ( m_InterpolatedValues[ 1 ] < 0 ) || ( m_InterpolatedValues[ 2 ] < 0 ) || ( m_InterpolatedValues[ 3 ] < 0 ) )
    {
    //std::cout << "pix < 0 kill" << std::endl;
    return true;
    }
  
  for ( int vertexNumber = 0; vertexNumber < 4; vertexNumber++ )
    {
    if ( m_InterpolatedValues[ vertexNumber ] == 0 )
      {
      if ( this->CheckBorderCase( m_NextRowAdditions[ vertexNumber ], 
                                  m_NextColumnAdditions[ vertexNumber ], 
                                  m_NextSliceAdditions[ vertexNumber ] ) )
        {
        return true;
        }
      }  
    }  

  // If we survived all these tests, we're inside
  return false;
    
}


//
//
//
template< typename TPixel >
bool
TetrahedronInteriorConstIterator< TPixel >
::CheckBorderCase( double a, double b, double c )
{ 
  if ( a < 0 )
    {
    return true;  
    }
    
  if ( a == 0 )
    {
    if ( b < 0 )
      {
      return true;  
      }
    if ( b == 0 )
      {
      if ( c < 0 )
        {
        return true;  
        }
      }
    }
    
 return false; 
  
}


} // end namespace kvl

#endif
