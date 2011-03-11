/**
 * @file  Broadcaster.cpp
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


#include <iostream>
#include "Broadcaster.h"
#include "Listener.h"

using namespace std;

Broadcaster::Broadcaster ( string isLabel ) :
    msLabel( isLabel ),
    m_bBlockBroadcast( false )
{
  mlListeners.clear();
}

Broadcaster::~Broadcaster ()
{
  for ( size_t i = 0; i < mlListeners.size(); i++ )
  {
    mlListeners[i]->RemoveBroadcaster( this );
  }
  mlListeners.clear();
}

void
Broadcaster::AddListener ( Listener* const iListener )
{
  for ( size_t i = 0; i < mlListeners.size(); i++ )
  {
    if ( iListener == mlListeners[i] )
      return;
  }
  mlListeners.push_back( iListener );
  iListener->AddBroadcaster( this );
}

void
Broadcaster::RemoveListener ( Listener* const iListener )
{
  for ( size_t i = 0; i < mlListeners.size(); i++ )
  {
    if ( iListener == mlListeners[i] )
    {
      mlListeners.erase( mlListeners.begin() + i );
      // iListener->RemoveBroadcaster( this );
      return;
    }
  }
}

void
Broadcaster::SendBroadcast ( std::string const iMessage,
                             void* iData, void* sender ) const
{
  if (m_bBlockBroadcast)
    return;

  for ( size_t i = 0; i < mlListeners.size(); i++ )
  {
    // cerr << "\tSending to listener " << listener->GetLabel() << endl;
    if ( !mlListeners[i]->IsBlocked() )
      mlListeners[i]->ListenToMessage( iMessage, iData, ( sender == NULL ? ((void*)this) : sender ) );
  }
}

void Broadcaster::BlockBroadcast( bool block )
{
  m_bBlockBroadcast = block;
}

