#ifndef __FieldLineInitPathStrategy_h
#define __FieldLineInitPathStrategy_h

#include "InitializePath.h"

#include <vnl/vnl_matrix.h>

class FieldLineInitPathStrategy : public InitializePath
{
public:

  FieldLineInitPathStrategy();
  FieldLineInitPathStrategy( PoistatsModel* model );
  ~FieldLineInitPathStrategy();  

  /**
   * This does the work of finding the initial path.
   */
  void CalculateInitialPath();
  
protected:

  /**
   * From the first and last seed point, calculate a look at vector and return
   * it.
   */
  void GetLookAtVector( double *lookAtVector );
  
  /**
   * Normalizes the vector and returns it.
   */
  void NormalizeVector( double *v, const int size );
  
  /**
   * Returns the up vector for an angle and look at vector.
   */
  void GetUpVector( double *up, const double *lookAt, const double angle );
  
  /**
   * Get the cross producet of 3 x 3 vectors.
   */
  void GetCrossProduct( double *product, const double *v1, const double *v2 );
  
  /**
   * From these vectors return the rotation matrix.
   */
  vnl_matrix< double > GetRotationMatrix( const double *lookAt, const double *up, 
    const double *right, const double* position );

};

#endif 

