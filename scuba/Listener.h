/**
 * @file  Listener.h
 * @brief Mixin class that recieves messages from Broadcasters
 *
 * Along with Broadcaster, this implements a broadcaster/listener
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


#ifndef Listener_h
#define Listener_h

#include "string_fixed.h"

class Listener {

public:

  // Ctor takes a label; not really used for anything other than
  // debugging.
  Listener ( std::string isLabel );
  virtual ~Listener ();

  std::string const& GetLabel() const;

  // Called by Broadcaster to receive a message.
  void ListenToMessage ( std::string isMessage, void* iData );

  // Subclasses should override this to respond to messages.
  virtual void DoListenToMessage ( std::string isMessage, void* iData ) = 0;

protected:

  std::string msLabel;
};

#endif
