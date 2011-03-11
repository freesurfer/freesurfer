/**
 * @file  Listener.cpp
 * @brief Mix-in class for message receiving
 *
 * Simple mix-in class for use with the Broadcaster class so text
 * messages with a pointer data can be sent to a list of listeners.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:39 $
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
#include "Listener.h"
#include "Broadcaster.h"

using namespace std;

Listener::Listener ( string isLabel ) :
    msLabel( isLabel ), m_bBlockListen( false )
{
  mBroadcasters.clear();
}

Listener::~Listener ()
{
  for ( size_t i = 0; i < mBroadcasters.size(); i++ )
  {
    mBroadcasters[i]->RemoveListener( this );
  }
  mBroadcasters.clear();
}

void Listener::ListenToMessage ( string const iMessage, void* iData, void* sender )
{

  //cerr << "Listener " << msLabel << " got message " << iMessage << endl;

  this->DoListenToMessage( iMessage, iData, sender );
}

void Listener::DoListenToMessage ( string const iMessage, void* iData, void* sender )
{

//  this->DoListenToMessage( iMessage, iData );
}

void Listener::BlockListen( bool bBlock )
{
  m_bBlockListen = bBlock;
}

bool Listener::IsBlocked()
{
  return m_bBlockListen;
}

void Listener::AddBroadcaster ( Broadcaster* const iBroadcaster )
{
  for ( size_t i = 0; i < mBroadcasters.size(); i++ )
  {
    if ( iBroadcaster == mBroadcasters[i] )
      return;
  }
  mBroadcasters.push_back( iBroadcaster );
}

void Listener::RemoveBroadcaster ( Broadcaster* const iBroadcaster )
{
  for ( size_t i = 0; i < mBroadcasters.size(); i++ )
  {
    if ( iBroadcaster == mBroadcasters[i] )
    {
      mBroadcasters.erase( mBroadcasters.begin() + i );
      return;
    }
  }
}
