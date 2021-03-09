/**
 * @brief icon...loader...test!
 *
 */
/*
 * Original Author: Nick Schmansky
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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
