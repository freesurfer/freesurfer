/**
 * @file  OrientMRIEvents.h
 * @brief Event constants
 *
 * These are the events used in our program.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/09/13 20:58:20 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2002-2007,
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

#ifndef OrientMRIEvents_h
#define OrientMRIEvents_h

#include "vtkCommand.h"

class OrientMRIEvents {

 public:

  // Description:
  // Our event constants.
  enum { VolumeToRASTransformChanged = vtkCommand::UserEvent + 1,
	 UserTransformChanged = vtkCommand::UserEvent + 2 };

};

#endif
