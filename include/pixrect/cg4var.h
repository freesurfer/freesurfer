/* @(#)cg4var.h 1.8 89/03/14 SMI */

/*
 * Copyright 1986 by Sun Microsystems, Inc.
 */

#ifndef	cg4var_DEFINED
#define	cg4var_DEFINED

/*
 * cg4 -- color memory frame buffer
 */

/* FBIOSATTR device specific array indices */
#define FB_ATTR_CG4_SETOWNER_CMD	0	/* 1 indicates PID is valid */
#define	FB_ATTR_CG4_SETOWNER_PID	1	/* new owner of device */

#include <pixrect/memvar.h>

#define	CG4_NFBS	3	/* number of frame buffers in a cg4 */

/* description of single frame buffer */
struct cg4fb {
	short			group;		/* plane group implemented */
	short			depth;		/* depth, bits */
	struct mprp_data	mprp;		/* memory pixrect data */
};
	
/* pixrect private data */
struct cg4_data {
	struct mprp_data	mprp;		/* memory pixrect simulator */
	int			flags;		/* misc. flags */
#define	CG4_PRIMARY		1		/* primary pixrect */
#define	CG4_OVERLAY_CMAP	2		/* overlay has colormap */
	int			planes;		/* current group and mask */
	int			fd;		/* file descriptor */
	short			active;		/* active fb no. */
	struct cg4fb		fb[CG4_NFBS];	/* frame buffer info */
};

/* useful macros */
#define	cg4_d(pr)	((struct cg4_data *) ((pr)->pr_data))

/* pixrect ops vector */
extern struct pixrectops cg4_ops;

int	cg4_putcolormap();
int	cg4_putattributes();
int	cg4_ioctl();

#ifndef KERNEL
Pixrect	*cg4_make();
int	cg4_destroy();
Pixrect *cg4_region();
int	cg4_getcolormap();
int	cg4_getattributes();
#endif !KERNEL

#endif	!cg4var_DEFINED
