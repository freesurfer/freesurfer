#ifndef _SBL_GEOMETRY_H_
#define _SBL_GEOMETRY_H_
#include <sbl/core/File.h>
#include <sbl/math/Vector.h>
namespace sbl {


//-------------------------------------------
// POINT 2 CLASS 
//-------------------------------------------


/// The Point2 class represents a 2D Euclidean point.
class Point2 {
public:

    // note: this class uses default copy constructor and assignment operator

    /// create a point
    Point2( double xIn, double yIn ) : x( xIn ), y( yIn ) {}

    /// coordinates
    double x;
    double y;

    // operators
    inline Point2 operator+( const Point2 &p ) const { return Point2( x + p.x, y + p.y ); }
};


//-------------------------------------------
// POINT 3 CLASS 
//-------------------------------------------


/// The Point3 class represents a 3D Euclidean point.
class Point3 {
public:

    // note: this class uses default copy constructor and assignment operator

    /// create a point
    Point3( double xIn, double yIn, double zIn ) : x( xIn ), y( yIn ), z( zIn ) {}
    Point3() : x( 0 ), y( 0 ), z( 0 ) {}

    /// coordinates
    double x;
    double y;
    double z;

    // operators
    inline Point3 operator+( const Point3 &p ) const { return Point3( x + p.x, y + p.y, z + p.z ); }
};


//-------------------------------------------
// SEGMENT 2 CLASS 
//-------------------------------------------


/// The Segment2 class represents a 2D line between two 2D points
class Segment2 {
public:

    // note: this class uses default copy constructor and assignment operator

    /// create a segment between two points
    Segment2( Point2 &start, Point2 &end ) : m_start( start ), m_end( end ) {}

    /// the segment start point
    const Point2 &start() const { return m_start; }

    /// the segment end point
    const Point2 &end() const { return m_end; }

private:

    // the points
    Point2 m_start;
    Point2 m_end;
};


//-------------------------------------------
// SEGMENT 3 CLASS 
//-------------------------------------------


/// The Segment3 class represents a 3D line between two 3D points
class Segment3 {

    // note: this class uses default copy constructor and assignment operator

    /// create a segment between two points
    Segment3( Point3 &start, Point3 &end ) : m_start( start ), m_end( end ) {}

    /// the segment start point
    const Point3 &start() const { return m_start; }

    /// the segment end point
    const Point3 &end() const { return m_end; }

private:

    // the points
    Point3 m_start;
    Point3 m_end;
};


//-------------------------------------------
// MATRIX 3 CLASS 
//-------------------------------------------


/// The Matrix3 class represents a 3x3 matrix for geometric transformations.
class Matrix3 {
public:

    // note: this class uses default copy constructor and assignment operator

    /// create a zero-initialized matrix
    Matrix3();

    /// the matrix data
    double data[ 3 ][ 3 ];

    /// set diagonal elements
    inline void setDiag( double x, double y, double z ) { data[ 0 ][ 0 ] = x; data[ 1 ][ 1 ] = y; data[ 2 ][ 2 ] = z; }

    /// set off-diagonal elements to a given value
    void setOffDiag( double v );

    /// transform a point using the matrix
    Point3 operator*( const Point3 &p ) const;
};


//-------------------------------------------
// AFFINE TRANSFORM 3 CLASS
//-------------------------------------------


/// The AffineTransform3 class represents a 3D linear transformation.
class AffineTransform3 {
public:

    // note: this class uses default copy constructor and assignment operator

    /// create identity transform
    AffineTransform3() { a.data[ 0 ][ 0 ] = 1; a.data[ 1 ][ 1 ] = 1; a.data[ 2 ][ 2 ] = 1; }

    /// create transform from vector of parameters
    AffineTransform3( const VectorD &params ); 

    /// the translation offset
    Point3 b;

    /// the transformation matrix (rotating, scaling, etc.)
    Matrix3 a;

    /// set the translation offset
    inline void setTranslation( double x, double y, double z ) { b.x = x; b.y = y; b.z = z; }

    /// set the diagonal matrix elements
    inline void setDiag( double x, double y, double z ) { a.data[ 0 ][ 0 ] = x; a.data[ 1 ][ 1 ] = y; a.data[ 2 ][ 2 ] = z; }

    /// set the off-diagonal matrix elements
    inline void setOffDiag( double v ) { a.setOffDiag( v ); }

    /// set the matrix elements
    void setMatrix( double a00, double a01, double a02, double a10, double a11, double a12, double a20, double a21, double a22 );

    /// set rotation matrix
    void setRotation( double x, double y, double z );

    /// convert the transformation to a vector of parameter values
    VectorD toVector() const;

    /// apply the linear transformation to a point
    inline Point3 transform( const Point3 &x ) const { return b + a * x; }

    /// compute the inverse transformation
    AffineTransform3 inverse() const;
};


/// save 3D affine transformation parameters to text file
void saveTransform( File &file, AffineTransform3 &transform );


/// load 3D affine transformation parameters from text file
AffineTransform3 loadTransform( File &file );


} // end namespace sbl
#endif // _SBL_GEOMETRY_H_

