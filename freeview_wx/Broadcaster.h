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
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:35 $
 *    $Revision: 1.1 $
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

#include <vector>
#include "string_fixed.h"

class Listener;

class Broadcaster
{

public:

  // Initialize with a string label, e.g. class name, that can be used
  // when debugging.
  Broadcaster ( std::string isLabel );
  virtual ~Broadcaster ();

  // Add or remove a listener from the list of listeners.
  void AddListener ( Listener* const iListener );
  void RemoveListener ( Listener* const iListener );

  void BlockBroadcast( bool block );

  // Send a message to all Listeners. They will have their
  // ListenToMessage (and subclassable DoListenToMessage) function
  // called with the string and pointer.
  virtual void SendBroadcast ( std::string const iMessage,
                               void* iData, void* sender = NULL ) const;

protected:
  std::string msLabel;
  std::vector<Listener*> mlListeners;

  bool m_bBlockBroadcast;
};


#endif
