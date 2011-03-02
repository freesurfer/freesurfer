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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.2 $
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
