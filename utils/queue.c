/**
 * @file  queue.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 01:49:39 $
 *    $Revision: 1.4 $
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
   @(#)queue.c  1.3
   3/15/94
*/
/*------------------------------------------------------------------------
      File Name: queue.c

         Author: Bruce Fischl

        Created: Jan. 1994

    Description:

------------------------------------------------------------------------*/
/*------------------------------------------------------------------------
                              HEADERS
------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "queue.h"
#include "thread.h"

/*------------------------------------------------------------------------
                            CONSTANTS
------------------------------------------------------------------------*/


/*------------------------------------------------------------------------
                            STATIC DATA
------------------------------------------------------------------------*/


/*------------------------------------------------------------------------
                            STATIC PROTOTYPES
------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
                              FUNCTIONS
------------------------------------------------------------------------*/
/*------------------------------------------------------------------------
       Parameters:

      Description:

    Return Values:
        nothing.
------------------------------------------------------------------------*/
QUEUE *
Qalloc(int max_elts)
{
  QUEUE *q ;
  int i=max_elts;

  q = (QUEUE *)calloc(1, sizeof(QUEUE)) ;

  i = max_elts+1 ;
  return(q) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:

    Return Values:
        nothing.
------------------------------------------------------------------------*/
void
Qfree(QUEUE *q)
{
  QELT *qelt ;

  for (qelt = Qget(q, 0) ; qelt ; qelt = Qget(q, 0))
    free(qelt) ;

  free(q) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:

    Return Values:
        nothing.
------------------------------------------------------------------------*/
int
Qput(QUEUE *q, void *data)
{
  QELT *qelt ;

  qelt = (QELT *)calloc(1, sizeof(QELT)) ;
  if (!qelt)
    return(-1) ;

  qelt->data = data ;

  if (!q->head)  /* empty list */
  {
    q->head = q->tail = qelt ;
    qelt->next = qelt->prev = NULL ;   /* just to be explicit */
  }
  else           /* link in at end of list */
  {
    qelt->next = NULL ;
    qelt->prev = q->tail ;
    q->tail->next = qelt ;
    q->tail = qelt ;
  }

  q->nelts++ ;
  return(0) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:

    Return Values:
        nothing.
------------------------------------------------------------------------*/
void *
Qget(QUEUE *q, int mode)
{
  QELT *qelt ;
  void *data ;

  if (!q->head)
  {
    if (mode == Q_DONT_WAIT)
      return(NULL) ;
    else
      ThreadSuspend(TID_SELF, 0) ;

    if (!q->head)   /* somebody woke me up and nothing was there */
      return(NULL) ;
  }

  qelt = q->head ;
  q->head = qelt->next ;
  if (!q->head)
    q->tail = NULL ;        /* empty list */
  else
    q->head->prev = NULL ;    /* head of list now */

  q->nelts-- ;
  data = qelt->data ;
  free(qelt) ;

  return(data) ;
}
