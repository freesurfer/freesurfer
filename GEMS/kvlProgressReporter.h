/**
 * @file  kvlProgressReporter.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:40 $
 *    $Revision: 1.3 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */
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
