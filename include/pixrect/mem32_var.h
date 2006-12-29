/* @(#)mem32_var.h 1.2 90/07/27 SMI */

/* Copyright 1990 Sun Microsystems, Inc. */

/* pixrect data for 8 bit emulation in 24/32 bit memory pixrect */

#ifndef  mem32var_DEFINED
#define  mem32var_DEFINED

#include <pixrect/pixrect.h>
#include <pixrect/memvar.h> /* mprp */
#include <sunwindow/cms.h> /* colormapseg */

#define         MEM32_8BIT_CMAPSIZE     256

struct mprp32_data
{
  struct mprp_data    mprp;
  int                 plane_group;
  int                 fd;
  int                 windowfd;
  struct colormapseg  cms;
};

#define  mprp32_d(pr)  ((struct mprp32_data *) (pr)->pr_data)

#ifndef KERNEL
extern
  struct pixrectops   mem32_ops;
int                 mem32_rop();
int                 mem32_getcolormap();
#endif /* !KERNEL */

#endif     /* !mem32var_DEFINED */
