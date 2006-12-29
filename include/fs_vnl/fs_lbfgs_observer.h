#ifndef FS_LBFGS_OBSERVER_H_
#define FS_LBFGS_OBSERVER_H_

#include <vnl/vnl_vector.h>

class fs_lbfgs_observer
{

public:

  fs_lbfgs_observer();
  ~fs_lbfgs_observer();

  int getNumberOfOptimalUpdates();

  void setStepFunction
  (
    void ( *stepFunction )
    ( int itno, float sse, void *parms, float *p ),
    void *params
  );

  void setUserCallbackFunction
  (
    void (*userCallbackFunction)(float [])
  );

  void update( double bestF, vnl_vector< double >* bestX );

private:

  int mNumberOfOptimalUpdates;

  void ( *mStepFunction )( int itno, float sse, void *parms, float *p );
  void *mStepFunctionParms;

  void ( *mUserCallbackFunction )( float []) ;

  void copyVnlToFloat( const vnl_vector<double>* input, float* output,
                       const int n);

  const bool hasStepFunction();

  const bool hasUserCallbackFunction();

};

#endif /*FS_LBFGS_OBSERVER_H_*/
