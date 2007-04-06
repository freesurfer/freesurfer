/**
 * @file  Listener.h
 * @brief Mix-in class for message receiving
 *
 * Simple mix-in class for use with the Broadcaster class so text
 * messages with a pointer data can be sent to a list of listeners.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/04/06 22:23:04 $
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


#ifndef Listener_h
#define Listener_h

#include "string_fixed.h"

class Listener {

public:

  // Initialize with a string label, e.g. class name, that can be used
  // when debugging.
  Listener ( std::string isLabel );
  virtual ~Listener ();

  // Called when a Broadcaster to which this object has ben added
  // through AddListener sends a message. This default is to call
  // DoListenToMessage.
  void ListenToMessage ( std::string const iMessage, void* const iData );

  // Subclasses should overrride this to respond to individual
  // messages.
  virtual void DoListenToMessage ( std::string const iMessage, void* const iData );

  std::string const GetLabel() const {
    return msLabel;
  }

protected:
  std::string msLabel;

};

#endif
