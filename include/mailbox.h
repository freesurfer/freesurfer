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
 *    $Date: 2011/03/02 00:04:09 $
 *    $Revision: 1.3 $
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
