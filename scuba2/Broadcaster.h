/**
 * @file  Broadcaster.h
 * @brief Mix-in class for message sending
 *
 * Simple mix-in class for use with the Listener class so text
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


#ifndef Broadcaster_h
#define Broadcaster_h

#include <list>
#include "string_fixed.h"
#include "Listener.h"

class Broadcaster {

public:

  // Initialize with a string label, e.g. class name, that can be used
  // when debugging.
  Broadcaster ( std::string isLabel );
  virtual ~Broadcaster ();

  // Add or remove a listener from the list of listeners. 
  void AddListener ( Listener* const iListener );
  void RemoveListener ( Listener* const iListener );

  // Send a message to all Listeners. They will have their
  // ListenToMessage (and subclassable DoListenToMessage) function
  // called with the string and pointer.
  virtual void SendBroadcast ( std::string const iMessage,
			       void* const iData=NULL ) const;

protected:
  std::string msLabel;
  std::list<Listener*> mlListeners;
};


#endif
