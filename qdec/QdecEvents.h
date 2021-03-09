/**
 * @brief Event constants
 *
 * These are the events used in our program.
 */
/*
 * Original Author: Kevin Teich
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

#ifndef QdecEvents_h
#define QdecEvents_h

#include "vtkCommand.h"

class QdecEvents {

 public:

  // Description:
  // Our event constants.
  enum { UserSelectedVertex = vtkCommand::UserEvent + 1 };

};

#endif
