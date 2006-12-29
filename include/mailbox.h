/**
 * @file  mailbox.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:08:59 $
 *    $Revision: 1.2 $
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


/*
   @(#)mailbox.h  1.2
   8/10/95
*/
/*------------------------------------------------------------------------
       File Name:  mailbox.h

         Author:  Bruce Fischl

        Created:  Jan. 1994

    Description:

------------------------------------------------------------------------*/
#ifndef MAILBOX_H
#define MAILBOX_H

#include <time.h>
#include <stdarg.h>
#include <sys/timeb.h>

typedef struct
{
  int    iMachineId ;      /* machine ID of sender */
  int    iTid ;            /* thread ID of sender */
  struct timeb timePosted ;      /* time message was posted */
  long   lMsgId ;          /* 4-character message ID */
  int    iUsers ;          /* # of users of this message */
  int    iLen ;            /* length of the data (excluding header) */
  void   *pvData ;
}
MSG ;

int  MBinit(int size) ;
int  MBprintf(long lMsgId, char *strFmt, ...) ;
int  MBsend(long iMsgId, int iLen, void *pvData) ;
int  MBreceive(long iMsgId, void (*func)(MSG *pmsg, void *parm),
               int iTid, void *parm) ;
int  MBinvoke(MSG *msg) ;


#endif
