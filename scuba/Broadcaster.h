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
