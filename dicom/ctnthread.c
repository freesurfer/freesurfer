/*
          Copyright (C) 1995 - 1996, RSNA and Washington University

          The software and supporting documentation for the Radiological
          Society of North America (RSNA) 1993 - 1996 Digital Imaging and
          Communications in Medicine (DICOM) Demonstration were developed
          at the
                  Electronic Radiology Laboratory
                  Mallinckrodt Institute of Radiology
                  Washington University School of Medicine
                  510 S. Kingshighway Blvd.
                  St. Louis, MO 63110
          as part of the 1993 - 1996 DICOM Central Test Node project for, and
          under contract with, the Radiological Society of North America.

          THIS SOFTWARE IS MADE AVAILABLE, AS IS, AND NEITHER RSNA NOR
          WASHINGTON UNIVERSITY MAKE ANY WARRANTY ABOUT THE SOFTWARE, ITS
          PERFORMANCE, ITS MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR
          USE, FREEDOM FROM ANY COMPUTER DISEASES OR ITS CONFORMITY TO ANY
          SPECIFICATION. THE ENTIRE RISK AS TO QUALITY AND PERFORMANCE OF
          THE SOFTWARE IS WITH THE USER.

          Copyright of the software and supporting documentation is
          jointly owned by RSNA and Washington University, and free access
          is hereby granted as a license to use this software, copy this
          software and prepare derivative works based upon this software.
          However, any distribution of this software source code or
          supporting documentation or derivative works (source code and
          supporting documentation) must include the three paragraphs of
          the copyright notice.
*/
/*
+-+-+-+-+-+-+-+-+-
*/
/*
**       Electronic Radiology Laboratory
**     Mallinckrodt Institute of Radiology
**  Washington University School of Medicine
**
** Module Name(s):
** Author, Date: Steve Moore, 30-Jun-96
** Intent:  Provide common abstractions needed for operations
**   in a multi-threaded environment.
** Last Update:  $Author: nicks $, $Date: 2006/12/29 02:08:57 $
** Source File:  $RCSfile: ctnthread.c,v $
** Revision:  $Revision: 1.5 $
** Status:  $State: Exp $
*/

#include <stdio.h>
#ifdef GCCSUNOS
#include <sys/types.h>
#endif
#include "dicom.h"
#include "ctnthread.h"

char *THR_Message(CONDITION cond);

#ifdef CTN_USE_THREADS
#ifdef _MSC_VER
#include <windows.h>
HANDLE hMutex[FAC_MAXIMUM];
#else
#include <synch.h>
static mutex_t mutex[FAC_MAXIMUM];
#endif
#endif

#ifdef CTN_USE_THREADS
static CTNBOOLEAN initialized = FALSE;
#endif

#ifdef _MSC_VER
static void
mutexName(int i, char *name) {
  sprintf(name, "CTN-MUTEX-%d", i);
}
#endif

CONDITION
THR_Init() {
#ifdef CTN_USE_THREADS
#ifdef _MSC_VER
  int i;

  if (initialized)
    return THR_NORMAL;

  for (i = 0; i < (int) DIM_OF(hMutex); i++) {
    char name[32];
    mutexName(i, name);
    hMutex[i] = CreateMutex(NULL, FALSE, name);
  }
  initialized = TRUE;
  return THR_NORMAL;
#else
  int cond;
  int i;

  if (initialized)
    return;

  for (i = 0; i < (int) DIM_OF(mutex); i++) {
    cond = mutex_init(&mutex[i], USYNC_THREAD, NULL);
    if (cond != 0) {
      fprintf(stderr, "Fatal error in THR_Init; could not initialize mutex\n");
      exit(1);
    }
  }
  initialized = TRUE;
  return THR_NORMAL;
#endif
#else
  return THR_NORMAL;
#endif
}

CONDITION
THR_Shutdown() {
#ifdef CTN_USE_THREADS
#ifdef _MSC_VER
  int i;

  if (!initialized) {
    fprintf(stderr, "Threads not initialized in call to THR_Shutdown\n");
    return THR_NOTINITIALIZED;
  }
  for (i = 0; i < (int) DIM_OF(hMutex); i++)
    CloseHandle(hMutex[i]);

  return THR_NORMAL;

#else
  int cond;
  int i;

  if (!initialized) {
    fprintf(stderr, "Threads not initialized in call to THR_Shutdown\n");
    return THR_NOTINITIALIZED;
  }
  for (i = 0; i < (int) DIM_OF(mutex); i++) {
    cond = mutex_destroy(&mutex[i]);
    if (cond != 0) {
      fprintf(stderr, "Failed on call to mutex_destroy in THR_Shutdown\n");
      return THR_GENERICFAILURE;
    }
  }
  return THR_NORMAL;
#endif
#else
  return THR_NORMAL;
#endif
}

CONDITION
THR_ObtainMutex(int fac) {
#ifdef CTN_USE_THREADS
#ifdef _MSC_VER
  char name[32];
  HANDLE hTmpMutex;
  mutexName(fac, name);
  hTmpMutex = OpenMutex(MUTEX_ALL_ACCESS, FALSE, name);
  WaitForSingleObject(hMutex[fac], INFINITE);
  /* From JCS, close temp handle to get eliminate resource leak. */
  CloseHandle(hTmpMutex);

  return THR_NORMAL;
#else
  int cond;

  if (!initialized) {
    fprintf(stderr,
            "Threads not initialized in call to THR_ObtainMutex: exiting\n");
    exit(1);
  }
  cond = mutex_lock(&mutex[fac]);
  if (cond != 0) {
    fprintf(stderr, "Failed on call to mutex_lock in THR_ObtainMutex\n");
    return THR_GENERICFAILURE;
  }
  return THR_NORMAL;
#endif
#else
  return THR_NORMAL;
#endif
}

CONDITION
THR_ReleaseMutex(int fac) {
#ifdef CTN_USE_THREADS
#ifdef _MSC_VER
  ReleaseMutex(hMutex[fac]);
  return THR_NORMAL;
#else
  int cond;

  if (!initialized) {
    fprintf(stderr, "Threads not initialized in call to THR_ReleaseMutex\n");
    return THR_NOTINITIALIZED;
  }
  cond = mutex_unlock(&mutex[fac]);
  if (cond != 0) {
    fprintf(stderr, "Failed on call to mutex_unlock in THR_ReleaseMutex\n");
    return THR_GENERICFAILURE;
  }
  return THR_NORMAL;
#endif
#else
  return THR_NORMAL;
#endif
}
