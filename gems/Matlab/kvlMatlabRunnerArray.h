#ifndef __kvlMatlabRunnerArray_h
#define __kvlMatlabRunnerArray_h

#include "kvlMatlabRunner.h"


namespace kvl
{


/**
  * Object to hold an array of Matlab runners.
  *
  * This class is a so-called Singleton, i.e. it is intelligent enough to 
  * ensure that only one single instance of it can be created. This is
  * useful in our case, as otherwise a new instance would be created and
  * populated with all the runners every time the Mex function is called.
  *
  */


class MatlabRunnerArray : public itk::Object
{
public:
  /** Smart pointer typedef support. */
  typedef MatlabRunnerArray         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(MatlabRunnerArray, itk::Object);

  /** Intercept calls to create new instances to ensure that we have a Singleton */
  static Pointer New();
  
  /** Return the singleton instance. */
  static Pointer GetInstance()
    { return New(); }

  /** */
  bool Run( const std::string& runnerName, 
            int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[] );
  
protected:
  MatlabRunnerArray();
  virtual ~MatlabRunnerArray() {};
  
  MatlabRunnerArray(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

  std::vector< MatlabRunner::Pointer >    m_Array;
  static Pointer  m_Instance;

};

} // end namespace kvl

#endif


