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
} MSG ;

int  MBinit(int size) ;
int  MBprintf(long lMsgId, char *strFmt, ...) ;
int  MBsend(long iMsgId, int iLen, void *pvData) ;
int  MBreceive(long iMsgId, void (*func)(MSG *pmsg, void *parm), 
                      int iTid, void *parm) ;
int  MBinvoke(MSG *msg) ;


#endif
