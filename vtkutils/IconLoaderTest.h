/**
 * @file  IconLoaderTest.h
 * @brief icon...loader...test!
 *
 */
/*
 * Original Author: Nick Schmansky
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/06/17 21:25:42 $
 *    $Revision: 1.1.2.1 $
 *
 * Copyright (C) 2008,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#ifndef __IconLoaderTest_h
#define __IconLoaderTest_h

#include "vtkKWApplication.h"
#include "vtkSmartPointer.h"

class vtkKWDialog;
class vtkKWTopLevel;

class IconLoaderTest : public vtkKWApplication {

 public:

  static IconLoaderTest* New ();

 protected:

  IconLoaderTest ();
  ~IconLoaderTest ();

};

#endif
