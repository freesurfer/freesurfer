#ifndef __kvlMatlabRunner_h
#define __kvlMatlabRunner_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "mex.h" 


namespace kvl
{


/**
  * Base class that defines the interface that
  * our Mex function can call to get work done
  *
  */


class MatlabRunner : public itk::Object
{
public:
  /** Smart pointer typedef support. */
  typedef MatlabRunner         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( MatlabRunner, itk::Object );

#if 0
  /** */
  virtual std::string GetName() const = 0;
#endif  
  
  /** */
  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {}
  
protected:
  MatlabRunner() {};
  virtual ~MatlabRunner() {};


  MatlabRunner(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl

#endif


