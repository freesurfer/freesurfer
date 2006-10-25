#include <iostream>

#include "fs_vnl/fs_lbfgs_observer.h"

fs_lbfgs_observer::fs_lbfgs_observer() {

  mNumberOfOptimalUpdates = 0;
  
  mStepFunction = NULL;
  mStepFunctionParms = NULL;
  mUserCallbackFunction = NULL;
  
}

fs_lbfgs_observer::~fs_lbfgs_observer() {
  
}

void fs_lbfgs_observer::update( double bestF, vnl_vector< double >* bestX ) {
    
  if( hasStepFunction() || hasUserCallbackFunction() ) {

    const int n = bestX->size();
    
      // legacy one indexing
    float currentX[n+1];
    copyVnlToFloat( bestX, currentX, n );
  
    if( hasStepFunction() ) {  
      (*mStepFunction)( mNumberOfOptimalUpdates, static_cast< float >( bestF ),
                        mStepFunctionParms, currentX);  
    }
    
    if( hasUserCallbackFunction() ) {  
      ( *mUserCallbackFunction )( currentX );
    }
    
  }

  mNumberOfOptimalUpdates++;
  
}

const bool fs_lbfgs_observer::hasStepFunction() {
  return ( mStepFunction != NULL );
}
    
const bool fs_lbfgs_observer::hasUserCallbackFunction() {
  return ( mUserCallbackFunction != NULL );
}

int fs_lbfgs_observer::getNumberOfOptimalUpdates() {
  return mNumberOfOptimalUpdates;
}

void fs_lbfgs_observer::setStepFunction
( void ( *stepFunction )
  ( int itno, float sse, void *parms, float *p ), void *parms)
{
  mStepFunction = stepFunction;
  mStepFunctionParms = parms;
}


void fs_lbfgs_observer::setUserCallbackFunction
  ( void (*userCallbackFunction)(float []) )
{
  mUserCallbackFunction = userCallbackFunction;
}

void fs_lbfgs_observer::copyVnlToFloat( const vnl_vector<double>* input, 
  float* output, const int n)
{
  for(int i=0; i<n; i++) {
    // legacy one indexing
    output[ i+1 ] = static_cast<float>( ( *input )( i ) );
  }
}

