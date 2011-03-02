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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:36 $
 *    $Revision: 1.7 $
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
