/**
 * @file  queue.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:00 $
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
   @(#)queue.h  1.4
   3/31/94
*/
/*------------------------------------------------------------------------
      File Name:  queue.h

         Author:  Bruce Fischl

        Created:  Jan. 1994

    Description:

------------------------------------------------------------------------*/


#ifndef QUEUE_H
#define QUEUE_H


typedef struct _qelt_
{
  struct _qelt_ *next ;
  struct _qelt_ *prev ;
  void          *data ;
}
QELT ;

typedef struct
{
  QELT  *head ;
  QELT  *tail ;
  int   nelts ;      /* # of elements in queue */
  int   mode ;
}
QUEUE ;

#define Q_WAIT_FOR_DATA   1
#define Q_DONT_WAIT       2

#ifdef ANSI
int   Qput(QUEUE *q, void *data) ;
void  *Qget(QUEUE *q, int mode) ;
QUEUE *Qalloc(int max_elts) ;
void  Qfree(QUEUE *q);
#else
int   Qput() ;
void  *Qget() ;
QUEUE *Qalloc() ;
void  Qfree();
#endif


#define Qempty(q)    (((q)->head == NULL))
#define Qfirst(q)    ((q)->head ? (q)->head->data : NULL)
#define Qnelts(q)    ((q)->nelts)

#endif



