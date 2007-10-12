/**
 * @file  Broadcaster.h
 * @brief Mixin class that sends messages to Listeners
 *
 * Along with Listener, this implements a broadcaster/listener
 * pattern. Classes can mixin the Broadcaster class and have classes
 * that have mixed in Listener added to their listener list, and then
 * send them all a string message with associated data.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/12 22:12:56 $
 *    $Revision: 1.6 $
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

  // Ctor takes a label; not really used for anything other than
  // debugging.
  Broadcaster ( std::string isLabel );
  virtual ~Broadcaster ();

  std::string const& GetLabel() const;

  // Add or remove a Listener from its list.
  void AddListener ( Listener& iListener );
  void RemoveListener ( Listener& iListener );

  // Send a message with optional data.
  virtual void SendBroadcast ( std::string isMessage, void* iData = NULL );

protected:

  std::string msLabel;
  std::list<Listener*> mlListeners;
};


#endif
