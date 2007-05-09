#include "FieldLineInitPathStrategy.h"

#include <iostream>

FieldLineInitPathStrategy::FieldLineInitPathStrategy() {
}

FieldLineInitPathStrategy::FieldLineInitPathStrategy( PoistatsModel* model ) { 
}

FieldLineInitPathStrategy::~FieldLineInitPathStrategy() { 
}  

void
FieldLineInitPathStrategy::CalculateInitialPath() {

  // this will be a vector with 3 columns
  std::vector< double* > path;
  
  // get the look at vector
  double lookAt[3];
  this->GetLookAtVector( lookAt );

  // TODO: we eventually want to rotate the up vector
  // use the look at as the normal to a plane that intersects a sphere to get
  // an up vector--a point on the unit circle
  const double upAngle = 0;
  double up[3];
  this->GetUpVector( up, lookAt, upAngle );
      
  // lets create circle with 5 points for the time being at z=0
  const int nPoints = 10;
  const double step = PI / static_cast< double >( nPoints );

  // TODO: remove this
  const int startLabel = ( *m_SeedValues )[ 0 ];
  int startPoint[3];
  this->GetRandomSeedPoint( startPoint, startLabel );

  // the position will the be translation component
  // TODO: what we really want is this to be the point along the spline
  const double position[] = { startPoint[0], startPoint[1], startPoint[2] };
  
  for( int i=0; i<nPoints; i++ ) {

    // add point to be saved
    double* point = new double[ 3 ];
    path.push_back( point );
    
    const double angle = static_cast< double >( i ) * step;
        
    // now that we have the up vector, we can calculate the right vector by
    // taking the cross product of the up and the look at
    double right[3];
    this->GetCrossProduct( right, lookAt, up );
        
    // create the rotation matrix
    vnl_matrix< double > rotation = this->GetRotationMatrix( lookAt, up, right, 
      position );
      
    // amplitude of the sine function
    const double amplitude = 5.0;
      
    // multiple point by the rotation
    vnl_matrix< double > rawPoint( 4, 1 );
    rawPoint( 0, 0 ) = i;
    rawPoint( 1, 0 ) = amplitude * sin( angle );
    rawPoint( 2, 0 ) = 0;
    rawPoint( 3, 0 ) = 1;
    
    vnl_matrix< double > transformedPoint = rotation * rawPoint;
    
    point[ 0 ] = transformedPoint( 0, 0 );
    point[ 1 ] = transformedPoint( 1, 0 );
    point[ 2 ] = transformedPoint( 2, 0 );

    std::cerr << "point: " << point[0] << ", " << point[1] << ", " << point[0] << std::endl;
    
  }

  // copy the path to the matrix output
  this->CopyPathToOutput( path );
  
}

void 
FieldLineInitPathStrategy::GetLookAtVector( double *lookAtVector ) {

  // get the vector pointing from the start point to the end
  const int startLabel = ( *m_SeedValues )[ 0 ];
  int startPoint[3];
  this->GetRandomSeedPoint( startPoint, startLabel );
  
  const int endLabel = ( *m_SeedValues )[ m_SeedValues->size()-1 ];
  int endPoint[3];
  this->GetRandomSeedPoint( endPoint, endLabel );
  
  // put the vector in our output vector
  for( int i=0; i<3; i++ ) {
    lookAtVector[ i ] = static_cast< double >( endPoint[ i ] - startPoint[ i ] );
  }
  
  // normalize the vector
  this->NormalizeVector( lookAtVector, 3 );
  
}

void 
FieldLineInitPathStrategy::NormalizeVector( double *v, const int size ) {
  
  // calculate the magnitude of the vector
  double magnitude = 0.0;
  for( int i=0; i<size; i++ ) {
    magnitude += v[ i ] * v[ i ];
  }
  
  magnitude = sqrt( magnitude );
  
  // normalize the vector based on the magnitude
  for( int i=0; i< size; i++ ) {
    v[ i ] /= magnitude;
  }
  
}

void 
FieldLineInitPathStrategy::GetUpVector( double *up, const double *lookAt, const double angle ) {
  
  // ideas from
  // http://en.wikipedia.org/wiki/Plane%E2%80%93sphere_intersection
  
  // lets treat the lookAt vector as the normal to a plane that cuts our sphere
  const double a = lookAt[ 0 ];
  const double b = lookAt[ 1 ];
  const double c = lookAt[ 2 ];
  
  // the equation that we're using for getting a unit vector is that of a point
  // on a cirlce that is the result of the intersection of a sphere and a plane
  
  double v1[] = { 0, c, -b };
  
  // normalize v1
  for( int i=0; i<3; i++ ) {
    v1[ i ] /= sqrt( 1 - a * a );
  }
  
  // v2 will be the cross product of v1 and the normal
  double v2[3];
  this->GetCrossProduct( v2, lookAt, v1 );
  
  // t and s are our parametric values for evaluating a point on a circle
  const double t = cos( angle );
  const double s = sin( angle );
  
  // up is a point on a circle
  for( int i=0; i<3; i++ ) {
    up[ i ] = s * v1[ i ] + t * v2[ i ];
  }
    
}

void 
FieldLineInitPathStrategy::GetCrossProduct( double *product, const double *v1, const double *v2 ) {

  product[0] = v1[ 1 ] * v2[ 2 ] - v1[ 2 ] * v2[ 1 ];
  product[1] = v1[ 2 ] * v2[ 0 ] - v1[ 0 ] * v2[ 2 ];
  product[2] = v1[ 0 ] * v2[ 1 ] - v1[ 1 ] * v2[ 0 ];
  
}

vnl_matrix< double >
FieldLineInitPathStrategy::GetRotationMatrix( const double *lookAt, const double *up, 
  const double *right, const double* position ) {
  
  vnl_matrix< double > rotation( 4, 4 );
  
  for( unsigned int cRow=0; cRow<rotation.rows(); cRow++ ) {
    
//    rotation( cRow, 0 ) = lookAt[ cRow ];
//    rotation( cRow, 1 ) = right[ cRow ];
//    rotation( cRow, 2 ) = up[ cRow ];

    // TODO: verify that this is right
    rotation( cRow, 0 ) = lookAt[ cRow ];
    rotation( cRow, 1 ) = up[ cRow ];
    rotation( cRow, 2 ) = right[ cRow ];

    rotation( cRow, 3 ) = position[ cRow ];
    
  }
  
  rotation( 3, 3 ) = 1.0;
  
  return rotation;
    
}
