/* @(#)cg12_var.h	1.12 of 5/9/91 SMI */

/* Copyright 1990 by Sun Microsystems, Inc. */

/* Sun Color Graphics board 12 (CG12) */

#ifndef cg12_var_DEFINED
#define	cg12_var_DEFINED

#include <sys/types.h>
#include <sys/ioctl.h>
#include <sun/fbio.h>
#include <pixrect/pixrect.h>
#include <sbusdev/cg12reg.h>
#include <pixrect/memvar.h>
#include <sunwindow/cms.h>		/* colormapseg */

#define	CG12_NFBS		6	/* number of frame buffers in a CG12 */

/* description of single CG12 frame buffer */

struct cg12fb
{
    short               group;		/* plane group implemented */
    short               depth;		/* depth, bits */
    struct mprp_data    mprp;		/* memory pixrect data */
};

/*  CG12  Pixrect private data */

struct cg12_data
{
    /* first make it look like a gp device */

    struct cg2fb       *cgpr_va;	/* backward source compatible */
    caddr_t             gp_shmem;	/* pointer to shared memory */
    int                 cgpr_fd;	/* primary flag */
    int                 cgpr_planes;	/* color bit plane mask reg */
    struct pr_pos       cgpr_offset;	/* pixrect offset */
    short               cg2_index;	/* cg2 board index */
    char                minordev;	/* true minor dev to stuff into GP */
    int                 gbufflag;	/* gbuffer flag */
    int                 ioctl_fd;	/* the fd to talk to the driver with */
    int                 ncmd;		/* length of cmdver array */
    u_char             *cmdver;		/* version #'s for each command */
    int                 flags;		/* misc options */
    int                 linebytes;	/* bytes per line (pixel mode) */
    int                 fbtype;		/* which cg is bound */

    /* then make it look like a memory device with multiple plane groups */

    struct mprp_data    mprp;		/* memory pixrect simulator */
    int                 cg12_flags;	/* misc. flags */
    int                 planes;		/* current group and mask */
    int                 fd;		/* file descriptor */
    short               active;		/* active fb no. */

    /* finally get some cg12 specific portions in */

    struct cg12_ctl_sp *ctl_sp;
    struct cg12fb       fb[CG12_NFBS];	/* frame buffer info */
    struct fb_wid_dbl_info	wid_dbl_info;	/* window id */
    int			windowfd;
};

/* cg12 specific constants */

/*
   There are 256 pixels of dither data.  The framebuffer is being addressed
   by an int ptr.  The framebuffer is in 4 bits/pixel mode.  Hence each 
   int points to 8 contiguous pixels.  Hence the offset from the base of 
   the framebuffer is  (num dith pixels + num solid pixels + num hollow
   pixels)/8.   See below.

   note, if the size of any of these sections is changed in the cg12
   microcode, then changes need to be propogated here too.  */

#define CG12_FB_SIZE   (CG12_WIDTH * CG12_HEIGHT)

/* amount to offset from beginning of offscreen memory for dither data */
#define CG12_DITHER_SIZE       256

/* amount to offset from beginning of offscreen memory for solid pattern data */
#define CG12_SOLID_PAT_SIZE    32

/* amount to offset from beginning of offscreen memory for hollow pattern data */
#define CG12_HOLLOW_PAT_SIZE   32

#define CG12_OFFSCREEN_START ((CG12_FB_SIZE + CG12_DITHER_SIZE + \
			       CG12_SOLID_PAT_SIZE + CG12_HOLLOW_PAT_SIZE) >> 3)

/* HACCESS value for 4 bits/pixel when loading offscreen VRAM */

#define CG12_OS_HACCESS 0x22

/* PLN_WR_MASK value for loading lower 4bits of each DPU offscreen VRAM */

#define CG12_OS_PLN_WR_MASK 0x000F0F0F

/* useful macros */

#define	cg12_d(pr)	((struct cg12_data *) ((pr)->pr_data))

#define CG12_PR_TO_MEM(src, mem)					\
    if (src && src->pr_ops != &mem_ops)					\
    {									\
	(void) cg12_set_state(src);					\
	mem.pr_ops      = &mem_ops;					\
	mem.pr_size     = src->pr_size;					\
	mem.pr_depth    = src->pr_depth;				\
	mem.pr_data     = (char *) &cg12_d(src)->mprp;			\
	src             = &mem;						\
    }

extern struct pixrectops cg12_ops;

int                 cg12_ioctl();
int                 cg12_putcolormap();
int                 cg12_putattributes();
int                 cg12_rop();

#ifndef KERNEL
extern int          gp1_rop();

int                 cg12_batchrop();
int                 cg12_destroy();
int                 cg12_get();
int                 cg12_getattributes();
int                 cg12_getcolormap();
Pixrect            *cg12_make();
int                 cg12_put();
Pixrect            *cg12_region();
int                 cg12_stencil();
int                 cg12_vector();

#endif	/* !KERNEL */

#endif	/* cg12var_DEFINED */
