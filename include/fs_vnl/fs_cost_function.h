#ifndef FS_COST_FUNCTION_H_
#define FS_COST_FUNCTION_H_

#define export // obsolete feature "export template" used in these header files
#include <vnl/vnl_cost_function.h>
#include <vnl/vnl_vector.h>
#undef export

class fs_cost_function : public vnl_cost_function
{

private:
  float (*mFunction)(float []);
  void (*mFunctionGradient)(float [], float []);

  void copyFromVNLToFloat
  ( float *floatVector,
    const vnl_vector< double > vnlVector, int numberOfParameters);


public:

  fs_cost_function( float (*function)(float []) );

  fs_cost_function( float (*function)(float []),
                    void (*functionGradient)(float [], float []),
                    int numberOfUnknowns );

  virtual ~fs_cost_function();

  void SetFunction( float (*function)(float []) );

  void SetFunctionGradient( void (*functionGradient)(float [], float []) );

  virtual double f(const vnl_vector<double> & x);

  virtual void gradf(const vnl_vector<double> & x,
                     vnl_vector<double>& gradient);

};

#endif /*FS_COST_FUNCTION_H_*/
