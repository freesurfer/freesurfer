/* @(#)tv1var.h 1.1 89/06/05 SMI */

/*
 * Copyright 1988 by Sun Microsystems, Inc.
 */

#ifndef tv1var_DEFINED
#define tv1var_DEFINED

#include <sys/types.h>
#include <sundev/tv1reg.h>
#include <pixrect/pixrect.h>
#include <pixrect/memvar.h>
#include <pixrect/cg4var.h>

#define TV1_NFBS 2
#define TV1_PRIMARY CG4_PRIMARY

struct tv1_data
{
  struct mprp_data mprp;
  int flags;
  int fd;
  short active;
  unsigned int planes;
  struct pr_pos offset;
  struct csr *tv1_csr;
  struct pr_size fbsize[TV1_NFBS];
  struct cg4fb    fb[TV1_NFBS];
  struct pixrect *emufb;
  struct pixrect *pr_video_enable;
  struct pixrect *pr_video;
};

struct tv1_enable_data
{
  struct pixrect *sub_pixrect;    /* pointer to the other pixrect */
  struct pr_pos  *offset;     /* pointer to offset, (in shared mem RO) */
  struct pr_pos  region_start;    /* start of region */
};


/*   DELETE THE FOLLOWING DURING CLEAN UP - ALL REFS SHOULD
     BE DIRECTLY VIA get_tv1_data -dw */
#define tv1_d(pr) ((struct tv1_data *)((pr)->pr_data))

#define tv1_enable_d(pr) ((struct tv1_enable_data *)((pr)->pr_data))

#ifndef KERNEL

Pixrect        *tv1_make ();
int  tv1_destroy ();
Pixrect        *tv1_region ();
int             tv1_getcolormap ();
int             tv1_getattributes ();

#ifndef ROUNDUP
#define ROUNDUP(val, gran) (((val) - 1 | (gran) - 1) + 1)
#endif

#endif !KERNEL

#endif tv1var_DEFINED
