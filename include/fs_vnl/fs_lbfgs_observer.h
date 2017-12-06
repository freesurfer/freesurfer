#ifndef FS_LBFGS_OBSERVER_H_
#define FS_LBFGS_OBSERVER_H_

#define export // obsolete feature "export template" used in these header files
#include <vnl/vnl_vector.h>
#undef export

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

  bool hasStepFunction();

  bool hasUserCallbackFunction();

};

#endif /*FS_LBFGS_OBSERVER_H_*/
