#include "FieldLineInitPathStrategy.h"

#include <iostream>

FieldLineInitPathStrategy::FieldLineInitPathStrategy( PoistatsModel* model ) {
  m_PoistatsModel = model;

  this->SetSeedVolume( m_PoistatsModel->GetSeedVolume() );
  this->SetSeedValues( m_PoistatsModel->GetSeedValues() );
  
  m_InitialPath = NULL;

  // this will be the seed value we use for getting random numbers
  std::srand( ( unsigned ) time( 0 ) );
  m_RandomTimeSeed = std::rand();
}

FieldLineInitPathStrategy::~FieldLineInitPathStrategy() { 
}  

void
FieldLineInitPathStrategy::CalculateInitialPath() {
  
  // this will be a vector with 3 columns
  std::vector< double* > path;

  // get the vector pointing from the start point to the end
  int startPoint[3];
  int endPoint[3];  
  this->GetNewStartAndEndPoints( startPoint, endPoint );

  // get the look at vector
  double lookAt[3];
  this->GetLookAtVector( lookAt, startPoint, endPoint );

  // TODO: we eventually want to rotate the up vector
  // use the look at as the normal to a plane that intersects a sphere to get
  // an up vector--a point on the unit circle
  // some random up angle
  const double upAngle = 2 * PI * m_PoistatsModel->GetRandomNumber();
  
  double up[3];
  this->GetUpVector( up, lookAt, upAngle );
  
  // number of times to sample the path
  // TODO: we probably only need the control points actually...
  const int nSamples = 25;
  
  // let the maximum amplitude be the radius from the model
  const double maxAmplitude = m_PoistatsModel->GetFieldLineRadius();
  
  // amplitude of the sine function
  const double amplitude = maxAmplitude * m_PoistatsModel->GetRandomNumber();

  // get the control points
  PoistatsModel::MatrixType initialPoints = this->GetInitialPoints( startPoint, endPoint );
  
  // rethread the path
  PoistatsModel::MatrixPointer rethreadedPath = 
    m_PoistatsModel->RethreadPath( &initialPoints, nSamples );    
  
  const double step = PI / static_cast< double >( nSamples );
  
  for( int i=0; i<nSamples; i++ ) {

    // the position will the be translation component
    const double position[] = { 
      ( *rethreadedPath )[ i ][ 0 ],
      ( *rethreadedPath )[ i ][ 1 ],
      ( *rethreadedPath )[ i ][ 2 ]
    };
  
    // add point to be saved
    double* point = new double[ 3 ];
    path.push_back( point );
    
    // now that we have the up vector, we can calculate the right vector by
    // taking the cross product of the up and the look at
    double right[3];
    this->GetCrossProduct( right, lookAt, up );
        
    // create the rotation matrix
    vnl_matrix< double > rotation = this->GetRotationMatrix( lookAt, up, right );
    
    // angle--position of the sine function
    const double angle = static_cast< double >( i ) * step;
        
    // multiple point by the rotation
    vnl_matrix< double > rawPoint( 3, 1 );
    rawPoint( 0, 0 ) = 0;
    rawPoint( 1, 0 ) = amplitude * sin( angle );
    rawPoint( 2, 0 ) = 0;
    
    vnl_matrix< double > transformedPoint = rotation * rawPoint;
    
    // translate the points
    for( int cRow=0; cRow<3; cRow++ ) {
      point[ cRow ] = transformedPoint( cRow, 0 ) + position[ cRow ];
    }

  }
  
  // we're missing the end point, so add it here
  double* point = new double[ 3 ];
  path.push_back( point );
  
  for( int nPoint=0; nPoint<3; nPoint++ ) {
    point[ nPoint ] = endPoint[ nPoint ];
  }

  // copy the path to the matrix output
  this->CopyPathToOutput( path );
  
  if( rethreadedPath != NULL ) {
    delete rethreadedPath;
    rethreadedPath = NULL;
  }
  
}

void 
FieldLineInitPathStrategy::GetLookAtVector( double *lookAtVector, 
  const int *startPoint, const int *endPoint ) {
  
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
  const double *right ) {
  
  vnl_matrix< double > rotation( 3, 3 );
  
  for( unsigned int cRow=0; cRow<rotation.rows(); cRow++ ) {
    
    // TODO: verify that this is right
    rotation( cRow, 0 ) = lookAt[ cRow ];
    rotation( cRow, 1 ) = up[ cRow ];
    rotation( cRow, 2 ) = right[ cRow ];

//    rotation( cRow, 3 ) = position[ cRow ];
    
  }
  
//  for( unsigned int cCol=0; cCol<rotation.cols(); cCol++ ) {
//    
//    // TODO: verify that this is right
//    rotation( 0, cCol ) = right[ cCol ];
//    rotation( 1, cCol ) = up[ cCol ];
//    rotation( 2, cCol ) = lookAt[ cCol ];
//    
//  }
  
  
  return rotation;
    
}

void 
FieldLineInitPathStrategy::GetNewStartAndEndPoints( int *startPoint, int *endPoint ) {
  
  std::vector< int > *seedValues = m_PoistatsModel->GetSeedValues();
  
  const int startLabel = ( *seedValues )[ 0 ];
  this->GetRandomSeedPoint( startPoint, startLabel );
  
  const int endLabel = ( *seedValues )[ seedValues->size()-1 ];
  this->GetRandomSeedPoint( endPoint, endLabel );
}

PoistatsModel::MatrixType 
FieldLineInitPathStrategy::GetInitialPoints( const int *startPoint, const int *endPoint ) {
  
  std::vector< int > *seedValues = m_PoistatsModel->GetSeedValues();
  
  PoistatsModel::MatrixType initialPoints( seedValues->size(), 3 );
  
  // set the start and end points
  for( unsigned int cCol=0; cCol<initialPoints.cols(); cCol++ ) {
    initialPoints( 0, cCol ) = startPoint[ cCol ];
    initialPoints( initialPoints.rows()-1, cCol ) = endPoint[ cCol ];
  }
  
  // set all the intermediate points if there are any
  for( unsigned int cRow=1; cRow<initialPoints.rows()-1; cRow++ ) {
    
    // get the intermediate point
    const int intermediateLabel = ( *seedValues )[ cRow ];
    int point[3];
    this->GetRandomSeedPoint( point, intermediateLabel );
    
    // set the point
    for( unsigned int cCol=0; cCol<initialPoints.cols(); cCol++ ) {
      initialPoints( cRow, cCol ) = point[ cCol ];
    }
    
  }
    
  return initialPoints;
  
}
