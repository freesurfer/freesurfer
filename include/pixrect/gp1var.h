/* @(#)gp1var.h 1.13 89/05/30 SMI */

/*
 * Copyright 1985, 1987 by Sun Microsystems, Inc.
 */

#ifndef gp1var_DEFINED
#define gp1var_DEFINED

/*
 * The same pixrect data structure is now used for both cg2
 * and gp1 pixrects.
 */
#include <sun/fbio.h>
#include <pixrect/cg2var.h>
#include <pixrect/cg9var.h>


struct gp1pr
{
  struct cg2fb *cgpr_va; /* backward source compatible */
  caddr_t  gp_shmem; /* pointer to shared memory */
  int   cgpr_fd; /* primary flag */
  int   cgpr_planes; /* color bit plane mask reg */
  struct pr_pos cgpr_offset; /* pixrect offset */
  short     cg2_index; /* cg2 board index */
  char  minordev; /* true minor dev to stuff into GP */
  int   gbufflag; /* gbuffer flag */
  int   ioctl_fd; /* the fd to talk to the driver with */
  int   ncmd;  /* length of cmdver array */
  u_char  *cmdver; /* version #'s for each command */
  int   flags;  /* misc options */
  int   linebytes; /* bytes per line (pixel mode) */
  int   fbtype;  /* which cg is bound */
  union {
    struct cg2pr  cg2pr;
    struct cg9_data cg9pr;
  } cgpr;
};

#define gp1_d(pr)  ((struct gp1pr *) (pr)->pr_data)
#define gp1_fbfrompr(pr) (gp1_d(pr)->cgpr.cg2pr.cgpr_va)
#define GP1_COPY_TO_CG2(gp_data, cg2_data) {                    \
 (cg2_data)->gp_shmem = (gp_data)->gp_shmem;             \
 (cg2_data)->cg2_index = (gp_data)->cg2_index;         \
 (cg2_data)->gbufflag = (gp_data)->gbufflag;             \
 (cg2_data)->minordev = (gp_data)->minordev;             \
 (cg2_data)->cmdver = (gp_data)->cmdver;                 \
 (cg2_data)->ncmd = (gp_data)->ncmd;                     \
 (cg2_data)->ioctl_fd = (gp_data)->ioctl_fd;             \
 (cg2_data)->flags = (gp_data)->flags;                   \
}

#define GP1_COPY_FROM_CG2(cg2_data, gp_data) {                  \
    (gp_data)->cgpr_fd = (cg2_data)->cgpr_fd;                   \
    (gp_data)->cgpr_planes = (cg2_data)->cgpr_planes;           \
    (gp_data)->cgpr_offset = (cg2_data)->cgpr_offset;           \
    (gp_data)->cgpr_va = (cg2_data)->cgpr_va;                   \
    (gp_data)->flags = (cg2_data)->flags;                       \
    (gp_data)->linebytes = (cg2_data)->linebytes;               \
}

#define GP1_COPY_TOTAL(srcdata, destdata) {                     \
    (destdata)->gp_shmem = (srcdata)->gp_shmem;                 \
    (destdata)->cgpr_fd = (srcdata)->cgpr_fd;                   \
    (destdata)->cgpr_va = (srcdata)->cgpr_va;                   \
    (destdata)->cgpr_planes = (srcdata)->cgpr_planes;           \
    (destdata)->cgpr_offset = (srcdata)->cgpr_offset;           \
    (destdata)->cg2_index = (srcdata)->cg2_index;               \
    (destdata)->minordev = (srcdata)->minordev;                 \
    (destdata)->gbufflag = (srcdata)->gbufflag;                 \
    (destdata)->ioctl_fd = (srcdata)->ioctl_fd;                 \
    (destdata)->ncmd = (srcdata)->ncmd;                         \
    (destdata)->cmdver = (srcdata)->cmdver;                     \
    (destdata)->flags = (srcdata)->flags;                       \
    (destdata)->linebytes = (srcdata)->linebytes;               \
}

struct gp1_version
{
  u_char majrel;
  u_char minrel;
  u_char serialnum;
  u_char flags;
};

#ifndef KERNEL
Pixrect *gp1_make();

extern struct pixrectops gp1_ops;

int gp1_rop();
int gp1_vector();
int gp1_destroy();
int gp1_put();
int gp1_get();
int gp1_putattributes();
int gp1_putcolormap();
int gp1_getattributes();
int gp1_getcolormap();
int gp1_ioctl();
Pixrect *gp1_region();
int gp1_batchrop();
int gp1_stencil();
int gp1_polypoint();

#endif !KERNEL

#define GP1IO_SATTR _IOW(G, 101, struct fbgattr)
#define GP1IO_SCMAP _IO(G, 102)
#define GP_SHMEMSIZE 5

#endif gp1var_DEFINED
