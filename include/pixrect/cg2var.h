/* @(#)cg2var.h 1.10 88/02/08 SMI */

/*
 * Copyright 1983, 1987 by Sun Microsystems, Inc.
 */

#ifndef	cg2var_DEFINED
#define	cg2var_DEFINED

/* FBIOSATTR device specific array indices */
#define	FB_ATTR_CG2_FLAGS	0	/* attribute flags word */
#define	FB_ATTR_CG2_FLAGS_SUN3	1	/* bit set for Sun-3 color board */
#define	FB_ATTR_CG2_FLAGS_PRFLAGS 2	/* bit set if PRFLAGS valid */

#define	FB_ATTR_CG2_BUFFERS	1	/* # of buffers */
#define	FB_ATTR_CG2_PRFLAGS	2	/* same as struct cg2pr flags */

struct cg2pr {
	struct cg2fb *cgpr_va;		/* mapped color board */
	caddr_t gp_shmem;		/* pointer to shared memory */
	int cgpr_fd;			/* primary flag */
	int cgpr_planes;		/* color bit plane mask reg */
	struct pr_pos cgpr_offset;	/* pixrect offset */
	short cg2_index;		/* cg2 board index */
	char minordev;			/* true minor dev to stuff into GP */
	int gbufflag;			/* gbuffer flag */
	int ioctl_fd;			/* the fd to talk to the driver with */
	int ncmd;			/* length of cmdver array */
	u_char *cmdver;			/* version #'s for each command */

	int flags;		/* misc options */
#define	CG2D_STDRES	1	/* standard (1152 x 900) resolution */
#define	CG2D_FASTREAD	2	/* has fast read feature */
#define	CG2D_ROPMODE	4	/* has aux. ropmode register */
#define	CG2D_32BIT	8	/* has 32 bit bus */
#define	CG2D_DBLBUF	16	/* has double buffering */
#define	CG2D_WINDMA	32	/* has window DMA */
#define	CG2D_ZOOM	256	/* has struct cg2_zoom */
#define	CG2D_NOZOOM	512	/* has struct cg2_nozoom */
#define	CG2D_ID		1024	/* has ID, extended status registers */
#define	CG2D_GP2	(1<<29)	/* GP is GP2 */
#define	CG2D_GB		(1<<30)	/* GB attached */
#define	CG2D_GP		(1<<31)	/* GP pixrect -- must be sign bit! */
	int linebytes;		/* bytes per line (pixel mode) */
};

#define cg2_d(pr) 	((struct cg2pr *)(pr)->pr_data)
#define cg2_fbfrompr(pr) (((struct cg2pr *)(pr)->pr_data)->cgpr_va)

#define cg2_ropword(cgd, plane, ax, ay)					\
		(cg2_ropwordaddr((cgd)->cgpr_va, (plane),		\
		 (cgd)->cgpr_offset.x+(ax),(cgd)->cgpr_offset.y+(ay)) )
#define cg2_pixel(cgd, ax, ay)						\
		(cg2_pixaddr((cgd)->cgpr_va,				\
		 (cgd)->cgpr_offset.x+(ax),(cgd)->cgpr_offset.y+(ay)) )
#define cg2_roppixel(cgd, ax, ay)					\
		(cg2_roppixaddr((cgd)->cgpr_va,				\
		 (cgd)->cgpr_offset.x+(ax),(cgd)->cgpr_offset.y+(ay)) )

#define cg2_prd_skew(cgd, ax)						\
		 (((cgd)->cgpr_offset.x+(ax)) & 15)

/* GP sync macro */
#define	GP1_PRD_SYNC(prd, error) _STMT( \
	if ((prd)->flags < 0 && gp1_sync((prd)->gp_shmem, (prd)->ioctl_fd)) \
		{ error; })

/* ops vector and functions */
extern	struct pixrectops cg2_ops;

int	cg2_rop();
int	cg2_putcolormap();
int	cg2_putattributes();

#ifndef KERNEL
int	cg2_stencil();
int	cg2_batchrop();
Pixrect *cg2_make();
int	cg2_destroy();
int	cg2_get();
int	cg2_put();
int	cg2_vector();
Pixrect *cg2_region();
int	cg2_getcolormap();
int	cg2_getattributes();
#endif !KERNEL

#endif	cg2var_DEFINED
