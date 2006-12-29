/* @(#)cg9var.h 1.10 90/03/22 SMI */

/* Copyright (c) 1989 by Sun Microsystems, Inc.  */

/*
 * CG9 Header File
 * Describe "Private Data" of CG9 Pixrect,  also external
 * references to pixrect operation
 */

#ifndef cg9var_DEFINED
#define cg9var_DEFINED

#include <sundev/cg9reg.h>
#include <pixrect/memvar.h>
#include <sunwindow/cms.h> /* colormapseg */

/* description of single CG9 frame buffer */

struct cg9fb
{
  short   group; /* plane group implemented */
  short   depth; /* depth, bits */
  struct mprp_data mprp; /* memory pixrect data */
};

#define CG_TYPE   3 /* for dev_specific stamp of fb.type */
#define CG_SIZE   4 /* for dev_specific stamp of fb size */
#define CG9_NFBS  3 /* number of frame buffers in a CG9 */

/*  CG9  Pixrect private data */

struct cg9_data
{
  struct mprp_data mprp; /* memory pixrect simulator */
  int   flags; /* misc. flags */
  int   planes; /* current group and mask */
  int   fd; /* file descriptor */
  short   active; /* active fb no. */
  struct cg9_control_sp *ctrl_sp;
  struct cg9fb  fb[CG9_NFBS]; /* frame buffer info */
  int   windowfd; /* if 8-bit indexed pw */
  struct colormapseg cms;  /* if 8-bit indexed pr */
  int   real_windowfd; /* for dblbuf hack */
};

/* useful macros */
#define cg9_d(pr) ((struct cg9_data *) ((pr)->pr_data))

#define CG9_PR_TO_MEM(src, mem)      \
    if (src && src->pr_ops != &mem_ops)     \
    {         \
 mem.pr_ops      = &mem_ops;     \
 mem.pr_size     = src->pr_size;     \
 mem.pr_depth    = src->pr_depth;    \
 mem.pr_data     = (char *) &cg9_d(src)->mprp;   \
 src             = &mem;      \
    }

#define CG9_PR_TO_MEM32(pr, mem32_pr, mem32_pr_data)   \
    if (pr && pr->pr_ops != &mem_ops)     \
    {         \
 mem32_pr.pr_ops = &mem32_ops;     \
 mem32_pr.pr_size = pr->pr_size;     \
 mem32_pr.pr_depth = pr->pr_depth;    \
 mem32_pr_data.mprp = cg9_d(pr)->mprp;    \
 mem32_pr_data.plane_group = cg9_d(pr)->fb[cg9_d(pr)->active].group;\
 mem32_pr_data.fd = cg9_d(pr)->fd;    \
 mem32_pr_data.windowfd = cg9_d(pr)->windowfd;   \
 mem32_pr_data.cms = cg9_d(pr)->cms;    \
 mem32_pr.pr_data = (char *) &mem32_pr_data;   \
 pr = &mem32_pr;       \
    }

extern struct pixrectops cg9_ops;
int  cg9_putcolormap();
int  cg9_putattributes();
int  cg9_rop();
int  cg9_ioctl();

#ifndef KERNEL
Pixrect  *cg9_make();
int  cg9_destroy();
Pixrect  *cg9_region();
int  cg9_getcolormap();
int  cg9_getattributes();
int  cg9_vector();
int  cg9_get();
int  cg9_put();

#endif !KERNEL

#endif cg9var_DEFINED
