/***************************************************************************
 *   Copyright (C) 2004 by Rudolph Pienaar / Christian Haselgrove          *
 *    Center for Morphometric Analysis       *
 *    Massachusetts General Hospital        *
 * Building 149, 13th St.         *
 *  Charlestown, MA 02129         *
 *     {ch|rudolph}@nmr.mgh.harvard.edu      *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
/// \file asynch.h
///
/// \brief Brief description
/// The asynchronous communications related object API.
///
/// \b DESCRIPTION
/// This defines the asynchronous communications parser and dispatching
/// layer.
///
/// \b HISTORY
/// 15 March 2005 - Initial consolidation from several other sources.
///
///

#ifndef __ASYNCH_H__
#define __ASYNCH_H__

#include "general.h"
#include "env.h"
#include "help.h"
#include "c_SSocket.h"

#ifdef __cplusplus
extern  "C" {
#endif

#include "mri.h"
#include "mrisurf.h"
#include "label.h"
#include "error.h"
#include "C_mpmProg.h"

#ifdef __cplusplus
}
#endif

#include <string>
using namespace std;

/// \fn void asynchEvent_poll(c_SSocket_UDP_receive* pCSocketUDPR, int maxPolls)
/// \brief The main socket polling function.
/// \param pCSocketUDPR  A socket on which to listen
/// \param maxPolls  The number of polls to wait. Each poll lasts a given timeout.
/// \return  (void)
string
asynchEvent_poll(
  c_SSocket_UDP_receive*  pCSocketUDPR,
  int    maxPolls
);


/// \fn void asynchEvent_process(s_env& st_env, string str_event)
/// \brief The main entry point to the asynchEvent process functions.
/// \param st_env   Core simulation environment
/// \param str_event  Control string parsed from socket
/// \return  (void)
void
asynchEvent_process(
  s_env&    st_env,
  string    str_event
);

/// \fn void asynchEvent_processSURFACE(s_env& st_env, string str_comms)
/// \brief Process socket-based events pertaining to the whole surface.
/// \param ast_env   Core simulation environment
/// \param astr_comms  Control string parsed from socket
/// \return  true|false  Status: success or fail.
bool
asynchEvent_processSURFACE(
  s_env&    st_env,
  string    astr_comms
);

/// \fn void asynchEvent_processLABEL(s_env& st_env, string str_comms)
/// \brief Process socket-based control of the LABEL object.
/// \param ast_env   Core simulation environment
/// \param astr_comms  Control string parsed from socket
/// \return  true|false  Status: success or fail.
bool
asynchEvent_processLABEL(
  s_env&    st_env,
  string    astr_comms
);

/// \fn void asynchEvent_processVERTEX(s_env& st_env, string str_comms)
/// \brief Process socket-based events pertaining to individual vertices.
/// \param ast_env   Core simulation environment
/// \param astr_comms  Control string parsed from socket
/// \return  true|false  Status: success or fail.
bool
asynchEvent_processVERTEX(
  s_env&    st_env,
  string    astr_comms
);

/// \fn void asynchEvent_processDWGHT(s_env& st_env, string str_comms)
/// \brief Process socket-based access to the Dweight structure
/// \param ast_env   Core simulation environment
/// \param astr_comms  Control string parsed from socket
/// \return  true|false  Status: success or fail.
bool
asynchEvent_processDWGHT(
  s_env&    ast_env,
  string    astr_comms
);

/// \fn void asynchEvent_processWGHT(s_env& st_env, string str_comms)
/// \brief Process socket-based access to the core weight structure
/// \param ast_env   Core simulation environment
/// \param astr_comms  Control string parsed from socket
/// \return  true|false  Status: success or fail.
bool
asynchEvent_processWGHT(
  s_env&    ast_env,
  string    astr_comms
);

/// \fn void asynchEvent_processENV(s_env& st_env, string str_comms)
/// \brief Process socket-based access to the core problem environment
/// \param ast_env   Core simulation environment
/// \param astr_comms  Control string parsed from socket
/// \return  true|false  Status: success or fail.
bool
asynchEvent_processENV(
  s_env&    ast_env,
  string    astr_comms
);

/// \fn void asynchEvent_processMPMPROG(s_env& st_env, string str_comms)
/// \brief Process socket-based access to the core problem environment
/// \param ast_env   Core simulation environment
/// \param astr_comms  Control string parsed from socket
/// \return  true|false  Status: success or fail.
bool
asynchEvent_processMPMPROG(
    s_env&      st_env,
    string      astr_comms
);

C_mpmProg_NOP*
pC_NOP_cast(
    C_mpmProg*                  pmpm,
    C_mpmProg_NOP*&             pC_mpmProg_NOP
);

C_mpmProg_autodijk*
pC_autodijk_cast(
    C_mpmProg*                  pmpm,
    C_mpmProg_autodijk*&        pC_mpmProg_autodijk
);

#endif //__ASYNCH_H__


