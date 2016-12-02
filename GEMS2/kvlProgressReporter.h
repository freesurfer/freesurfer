#ifndef __kvlProgressReporter_h
#define __kvlProgressReporter_h

#include "itkSimpleFilterWatcher.h"


/**
 *
 */
namespace kvl
{

class ProgressReporter : public itk::SimpleFilterWatcher
{
public:

 //
 ProgressReporter( itk::ProcessObject* o, const char *comment="" )
   : SimpleFilterWatcher( o, comment ) {}

 ProgressReporter( const itk::ProcessObject* o, const char *comment="" )
   : SimpleFilterWatcher( const_cast< itk::ProcessObject* >( o ), comment ) {}

protected:

  // Overload some stuff to only show what we need
  virtual void ShowProgress()
    {
    std::cout << "  " << this->GetProcess()->GetProgress() * 100 << "%" << std::endl;
    }

  virtual void StartFilter()
    {
    std::cout << this->GetComment() << ": " << std::endl;
    }

  virtual void EndFilter()
    {
    std::cout << "  100%" << std::endl;
    }

private:

};


} // end namespace kvl

#endif
