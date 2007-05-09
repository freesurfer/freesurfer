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
  
  // lets create circle with 5 points for the time being at z=0
  const int nPoints = 5;
  const double step = PI / nPoints;
  
  const double z = 0;
  for( int i=0; i<nPoints; i++ ) {

    // add point to be saved
    double* point = new double[ 3 ];
    path.push_back( point );
    
    const double angle = static_cast< double >( i ) * step;

    const double x = sin( angle );
    const double y = cos( angle );
    
    std::cerr << "point: " << x << ", " << y << ", " << z << std::endl;
    
    point[ 0 ] = x;
    point[ 1 ] = y;
    point[ 2 ] = z;
    
  }

  // copy the path to the matrix output
  this->CopyPathToOutput( path );
  
}
