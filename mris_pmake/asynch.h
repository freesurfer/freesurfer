/**
 * @brief The asynchronous communications related object API
 *
 * This defines the asynchronous communications parser and dispatching layer.
 */
/*
 * Original Author: Rudolph Pienaar / Christian Haselgrove
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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


#ifndef __ASYNCH_H__
#define __ASYNCH_H__

#include "general.h"
#include "env.h"
#include "help.h"
#include "c_SSocket.h"

#include "mri.h"
#include "mrisurf.h"
#include "label.h"
#include "error.h"
#include "C_mpmProg.h"

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

#if 0
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
#endif

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

C_mpmProg_pathFind*
pC_pathFind_cast(
    C_mpmProg*                  pmpm,
    C_mpmProg_pathFind*&        pC_mpmProg_pathFind
);

C_mpmProg_ROI*
pC_ROI_cast(
    C_mpmProg*                  pmpm,
    C_mpmProg_ROI*&             pC_mpmProg_ROI
);

C_mpmProg_externalMesh*
pC_externalMesh_cast(
    C_mpmProg*                  pmpm,
    C_mpmProg_externalMesh*&	pC_mpmProg_exernalMesh
);


#endif //__ASYNCH_H__


